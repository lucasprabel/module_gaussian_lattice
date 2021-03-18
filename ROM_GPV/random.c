#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <fcntl.h>
#include <unistd.h>
#include <x86intrin.h>

#include "common.h"
#include "random.h"

/*
	Returns an integer sampled from the uniform distribution over [0, n]
		using the uniform distribution over [0, 2^32 - 1] provided by random_bytes
*/
uint32_t uniform_int_distribution(uint32_t n)
	{
	uint32_t scaling = (UINT32_MAX) / (n + 1);
	uint32_t past = (n + 1) * scaling;
	uint32_t r_data[4];
	
	while(true)
		{
		random_bytes((uint8_t *) r_data);
		
		for(int i = 0 ; i < 4 ; ++i)
			{
			uint32_t r = r_data[i];
			
			if(r < past)
				{
				return r / scaling;
				}
			}
		}
	}

/*
  Returns an integer sampled from the uniform distribution over [0, q-1]
    using the uniform distribution over [0, 2^32 - 1] provided by random_bytes
*/
scalar uniform_mod_q(void)
  {
  uint32_t scaling = (UINT32_MAX) / PARAM_Q;
  uint32_t past = PARAM_Q * scaling;
  uint32_t r_data[4];
  
  while(true)
    {
    random_bytes((uint8_t *) r_data);
    
    for(int i = 0 ; i < 4 ; ++i)
      {
      uint32_t r = r_data[i];
      
      if(r < past)
        {
        return r / scaling;
        }
      }
    }
  }
	
#define AES_128_key_exp(k, rcon) aes_128_key_expansion(k, _mm_aeskeygenassist_si128(k, rcon))
#define DO_ENC_BLOCK(m,k) \
  do {\
        m = _mm_xor_si128       (m, k[ 0]); \
        m = _mm_aesenc_si128    (m, k[ 1]); \
        m = _mm_aesenc_si128    (m, k[ 2]); \
        m = _mm_aesenc_si128    (m, k[ 3]); \
        m = _mm_aesenc_si128    (m, k[ 4]); \
        m = _mm_aesenc_si128    (m, k[ 5]); \
        __auto_type m5 = m;\
        m = _mm_aesenc_si128    (m, k[ 6]); \
        m = _mm_aesenc_si128    (m, k[ 7]); \
        m = _mm_aesenc_si128    (m, k[ 8]); \
        m = _mm_aesenc_si128    (m, k[ 9]); \
        m = _mm_aesenclast_si128(m, k[10]);\
        m = _mm_xor_si128(m, m5);\
    } while(0)

//the expanded key
static __m128i key_schedule[11];
static uint64_t ctr = 0;


static __m128i aes_128_key_expansion(__m128i key, __m128i keygened) {
  keygened = _mm_shuffle_epi32(keygened, _MM_SHUFFLE(3,3,3,3));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  return _mm_xor_si128(key, keygened);
}

static void aes128_load_key(const int8_t * enc_key){
    key_schedule[0] = _mm_loadu_si128((const __m128i*) enc_key);
  key_schedule[1]  = AES_128_key_exp(key_schedule[0], 0x01);
  key_schedule[2]  = AES_128_key_exp(key_schedule[1], 0x02);
  key_schedule[3]  = AES_128_key_exp(key_schedule[2], 0x04);
  key_schedule[4]  = AES_128_key_exp(key_schedule[3], 0x08);
  key_schedule[5]  = AES_128_key_exp(key_schedule[4], 0x10);
  key_schedule[6]  = AES_128_key_exp(key_schedule[5], 0x20);
  key_schedule[7]  = AES_128_key_exp(key_schedule[6], 0x40);
  key_schedule[8]  = AES_128_key_exp(key_schedule[7], 0x80);
  key_schedule[9]  = AES_128_key_exp(key_schedule[8], 0x1B);
  key_schedule[10] = AES_128_key_exp(key_schedule[9], 0x36);
}

void random_bytes_init(void) {
    int8_t seed[16] = {0x01, 0x05, 0x07, 0x09, 0x0F, 0xFF};
    int randomData = open("/dev/urandom", O_RDONLY);
    read(randomData, seed, sizeof(int8_t)*16);
//    _rdseed64_step((unsigned long long * ) seed);
//    _rdseed64_step((unsigned long long * ) (seed + 8));
    aes128_load_key(seed);
    ctr = 0;
}

void random_bytes(uint8_t * restrict data) {
    __m128i m = _mm_setr_epi32(0, 1, 2, ctr++);
    DO_ENC_BLOCK(m, key_schedule);

    uint8_t * x = (uint8_t *) __builtin_assume_aligned(data, 16);
    _mm_storeu_si128((__m128i *) x, m);
}

static const double ln2 = 0.6931471805599453; // log(2);
static const double a = 5.713363152645423; //(4+3*sqrt(2))*log(2);
static const double b = 3.414213562373095; // 2+sqrt(2);
static const double c = -1.6734053240284923; // -(1+sqrt(2))*log(2);
static const double p = 0.9802581434685472; //sqrt(2)*log(2);
static const double A = 5.6005707569738075; //a*p;
static const double B = 3.3468106480569846; // b*p;
static const double D = 0.08578643762690495; // 1/(b*b);
static const double H = 0.0026106723602095233;// h*D/p;
int algorithmF(const double mu, const double sigma) {
  static uint32_t randomData[8];
  static unsigned randomDataIndex = 8;

  static uint8_t randomExpSecond[32];
  static unsigned randomExpSecondIndex = 8;

  const double _sigma_floor_ = (unsigned) sigma;  
  const unsigned sigma_floor = (unsigned) (((_sigma_floor_ == sigma) && (_sigma_floor_ != 0)) ? sigma - 1 : (_sigma_floor_));

  const uint32_t scaling = UINT32_MAX / (sigma_floor + 1);
  const uint32_t past = (sigma_floor + 1) * scaling;
  
  while(true) {
    do {
      if(randomDataIndex == 8) {
        randomDataIndex = 0;
        random_bytes((uint8_t *) randomData);
        random_bytes((uint8_t *) (randomData + 4));
      }

    } while(randomData[randomDataIndex++] >= past);

    
    // compute j and s
    const unsigned j = randomData[randomDataIndex - 1] / scaling;
    const int s = ((randomData[randomDataIndex - 1] & 0x01) == 0) ? -1 : 1;
        

    // Exp
    double expVariate[3];
    double expG[3];
    unsigned expVariateIndex = 0;
    
    // First try (high probability)
    uint8_t randomExpFirst[32];
    random_bytes(randomExpFirst);
    random_bytes(randomExpFirst + 16);

    const double U_1 = *((uint16_t *) randomExpFirst) / (UINT16_MAX*1.0);
    {
      const unsigned long long int I_i = *((uint64_t *) (randomExpFirst + 2));
      const int j_i = __builtin_ctzll(I_i);
      expG[expVariateIndex] = c + j_i*ln2;
    }   
    if (U_1 <= p) {
      expVariate[expVariateIndex] = expG[expVariateIndex] + A/(B - U_1);
      expVariateIndex++;
    }

    const double U_2 = *((uint16_t *) (randomExpFirst + 10)) / (UINT16_MAX*1.0);
    {
      const unsigned long long int I_i = *((uint64_t *) (randomExpFirst + 12));
      const int j_i = __builtin_ctzll(I_i);
      expG[expVariateIndex] = c + j_i*ln2;
    }
    
    if (U_2 <= p) {
      expVariate[expVariateIndex] = expG[expVariateIndex] + A/(B - U_2);
      expVariateIndex++;
    }
    
    const double U_3 = *((uint16_t *) (randomExpFirst + 20)) / (UINT16_MAX*1.0);
    {
      const unsigned long long int I_i = *((uint64_t *) (randomExpFirst + 22));
      const int j_i = __builtin_ctzll(I_i);
      expG[expVariateIndex] = c + j_i*ln2;
    }
    if (U_3 <= p) {
      expVariate[expVariateIndex] = expG[expVariateIndex] + A/(B - U_3);
      expVariateIndex++;
    }


    for(; expVariateIndex < 3;) {
      for(; randomExpSecondIndex < 8; randomExpSecondIndex++) {
        const double U1 = *((uint16_t *) (randomExpSecond + randomExpSecondIndex*2)) / (UINT16_MAX*1.0);
        const double U2 = *((uint16_t *) (randomExpSecond + randomExpSecondIndex*2 + 8)) / (UINT16_MAX*1.0);
        const double bU1 = b - U1;
        const double Y = a/bU1;
        const double L = (U2*H + D)*(bU1)*(bU1);
        const double Z = Y + c;

        const double LZ  = L - 1 + Z;
        const double ZZ  = Z*Z;
        const bool cond1 = LZ <= 0;
        const bool cond2 = LZ - ZZ/2 + ZZ*Z/6 <= 0;
        const bool cond3 = LZ - ZZ <= 0;


        if (cond1 || cond2 || (cond3 && (L <= exp(-Z)))) {
          expVariate[expVariateIndex] = expG[expVariateIndex] + Y;
          expVariateIndex++;
          break;
        }
      }
      if(randomExpSecondIndex == 8) {
        randomExpSecondIndex = 0;
        random_bytes(randomExpSecond);
        random_bytes(randomExpSecond + 16);
      } 
    }

    unsigned k = ((unsigned) (2*expVariate[0]));
    int i0 = ((unsigned) (sigma*k + s*mu)) + 1;
    double x0 = (i0 - (sigma*k + s*mu))/sigma;
    double x = x0 + j/sigma;

    unsigned k1 = ((unsigned) (2*expVariate[1]));
    double z = expVariate[2];

    bool cond1 = k1 >= k*(k-1);
    bool cond2 = x < 1;
    bool cond3a = x != 0;
    bool cond3b = k != 0;
    bool cond3c = s == 1;
    bool cond3 = cond3a || cond3b || cond3c;
    bool cond4 = z > 0.5*x*(2*k + x);

    if (cond1 && cond2 && cond3 && cond4) {
      return s*(i0 + j);
    }
  }
}

/*
	Generates a random salt r of size SALT_BYTES
*/
void salt(uint8_t *r)
	{
	for(int i = 0 ; i < NEEDED_AES_FOR_SALT ; ++i)
		{
		random_bytes(&r[16*i]);
		}
	}
