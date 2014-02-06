#ifndef __CONFIG_H__
#define __CONFIG_H__


#include <stdio.h>


#define _DEBUG_LEVEL_    1   //  0 to 10,
                             //  0 is no debug print info at all,
                             //  10 is full info

#if defined (__MIC__)
#define SIMDW      64
#elif defined (__AVX__)
#define SIMDW      32
#elif defined (__SSE__)
#define SIMDW      16
#else
#define SIMDW      64
#endif


#define ALIGNED_MALLOC(size)  _mm_malloc(size, __ALIGNLEN__)
#define ALIGNED_FREE(addr)    _mm_free(addr)


#if ( _DEBUG_LEVEL_ == -1 )
#define CINT_PRINTF( level, fmt, args... )        {}
#else
#define CINT_PRINTF( level, fmt, args... )                                          \
        do                                                                          \
        {                                                                           \
            if ( (unsigned)(level) <= _DEBUG_LEVEL_ )                               \
            {                                                                       \
                sprintf( basis->str_buf, "%s() line %d ", __FUNCTION__, __LINE__ ); \
                sprintf( basis->str_buf + strlen(basis->str_buf), fmt, ##args );    \
                fprintf( stdout, "%s", basis->str_buf );                            \
                fflush( stdout );                                                   \
            }                                                                       \
        } while ( 0 )
#endif


#if ( _DEBUG_LEVEL_ > 1 )
#define CINT_INFO( fmt, args... )                                            \
        do                                                                   \
        {                                                                    \
            sprintf( basis->str_buf, "**** CInt: ");                         \
            sprintf( basis->str_buf + strlen("**** CInt: "), fmt, ##args );  \
            fprintf( stdout, "%s", basis->str_buf );                         \
            fflush( stdout );                                                \
        } while ( 0 )
#else
#define CINT_INFO( fmt, args... )        {}
#endif


#endif /* __CONFIG_H__ */
