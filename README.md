# LwFEC - Lightweight Forward Error Correction library
LwFEC is a full Reed-Solomon FEC encoder and decoder C library dedicated for embedded and other RAM usage restricted systems.
## Memory usage
 The aim of this library is to be safe and deterministic, thus:
* No heap (malloc) is used
* No stack-allocated arrays are used

Block size limit and parity byte count limit can be changed freely according to application requirements. This allows the memory usage to stay at the required minimum.
There is an architectural limit of 256 block size limit due to use of GF(2^8).

## Example usage
The example below shows how to create generator polynomial, encode and decode a message:
```C
#include <stdint.h> //standard integer types header
#include "rs.h" //RS FEC library

#define N 16 //block (data + parity) size
#define T 5 //parity size

uint8_t generator[T + 1]; //generator polynomial buffer

void LwFEC(void)
{
    RsGeneratePolynomial(T, generator); //create generator polynomial
    uint8_t data[N] = {'H', 'e', 'l', 'l', 'o', ' ', 'w', 'o', 'r', 'l', 'd', 0, 0, 0, 0, 0, 0};
    RsEncode(data, N, generator, T); //encode message

    //then N bytes of data are sent through some channel
    //or stored on some disk
    //this may introduce errors that the decoder will try to fix

    uint8_t bytesFixed = 0; //store number of corrected bytes here
    if(RsDecode(data, N, generator, T, &bytesFixed)) //decode message
    {
        //message decoded succesfully
    }
    else
    {
        //too many errors in message
    }
}
```
## Important notes
The library does not provide any kind of data segmentation, padding etc. The coder/decoder does not store any configuration information, thus it can be used dynamically with any block size *N* and data size *K* (or parity size *T=N-K*), as long as the appropriate generator polynomial is provided.
## Size limit change
Size limits used for array preallocation are stored in *rs.h*:
```C
#define RS_MAX_BLOCK_SIZE 256 //maximum block (data + parity bytes) size
#define RS_MAX_REDUNDANCY_BYTES 32 //maximum parity bytes
```
## License
The code is based on [*Reed-Solomon codes for coders*](https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders).
The project is licensed under the GNU GPL v3 license (see [LICENSE](LICENSE)).

