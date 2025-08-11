/*
Copyright 2023 Piotr Wilkoń
This file is part of LwFEC.

LwFEC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

LwFEC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LwFEC.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include <stdbool.h>

#define RS_MAX_REDUNDANCY_BYTES 64 //maximum parity bytes


#define RS_BLOCK_SIZE 255 //natural full RS block size
#define RS_MAX_DATA_SIZE (RS_BLOCK_SIZE - RS_MAX_REDUNDANCY_BYTES)

/**
 * @brief Reed-Solomon module configuration and operation structure
*/
struct LwFecRS
{
    uint8_t generator[RS_MAX_REDUNDANCY_BYTES + 1]; //generator polynomial
    uint8_t T; //number of redundancy/parity bytes
    uint8_t fcr; //first consecutive root index


    union
    {
        uint8_t initBuffer[RS_MAX_REDUNDANCY_BYTES + 1];
        uint8_t encodeBuffer[RS_BLOCK_SIZE];
    };

    struct
    {
        uint8_t syndromes[RS_MAX_REDUNDANCY_BYTES + 1];
        uint8_t locator[RS_MAX_REDUNDANCY_BYTES + 1];
        uint8_t evaluator[RS_MAX_REDUNDANCY_BYTES + 1];
    } decodeBuffer;

    union
    {
        union
        {
            struct
            {
                uint8_t errLoc[RS_MAX_REDUNDANCY_BYTES + 1];
                uint8_t newLoc[RS_MAX_REDUNDANCY_BYTES + 1];
                uint8_t oldLoc[RS_MAX_REDUNDANCY_BYTES + 1];
                uint8_t tmpLoc[RS_MAX_REDUNDANCY_BYTES + 1];
            };
            struct
            {
                uint8_t B[RS_MAX_REDUNDANCY_BYTES + 1];
                uint8_t C[RS_MAX_REDUNDANCY_BYTES + 1];
                uint8_t T[RS_MAX_REDUNDANCY_BYTES + 1];
                uint8_t T2[RS_MAX_REDUNDANCY_BYTES + 1];
            };
        } locatorBuffer;
        struct
        {
            uint8_t locator[RS_MAX_REDUNDANCY_BYTES + 1];
            uint8_t errataEvaluator[2 * RS_MAX_REDUNDANCY_BYTES + 2];
            uint8_t errataPosition[RS_MAX_REDUNDANCY_BYTES];
        } fixerBuffer;
        
    };
    
};

/**
 * @brief Decode message using Reed-Solomon FEC
 * 
 * This function takes input buffer with K data bytes and T parity bytes
 * Then it moves parity bytes to the end and fills everything inbetween with zeros.
 * Next the in-place decoding is performed.
 * @param *rs RS coder/decoder instance
 * @param *data Input/output buffer. Must be of size N = 255
 * @param size Data size = K
 * @param *fixed Output number of bytes corrected
 * @return True on success, false on failure
 */
bool RsDecode(struct LwFecRS *rs, uint8_t *data, uint8_t size, uint8_t *fixed);

/**
 * @brief Encode message using Reed-Solomon FEC
 * @param *rs RS coder/decoder instance
 * @param *data Input/output buffer. Must be of size N = 255
 * @param size Data size = K
 */
void RsEncode(struct LwFecRS *rs, uint8_t *data, uint8_t size);

/**
 * @brief Initialize Reed-Solomon coder/decoder
 * 
 * This function calculates generator polynomial and stores required constants.
 * @param *rs RS coder/decoder instance to be filled
 * @param T Number of parity check bytes
 * @param fcr First consecutive root index
*/
void RsInit(struct LwFecRS *rs, uint8_t T, uint8_t fcr);
