/*
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

#define RS_MAX_BLOCK_SIZE 255 //maximum block (data + parity bytes) size
#define RS_MAX_REDUNDANCY_BYTES 64 //maximum parity bytes

#define RS_MAX_DATA_SIZE (RS_MAX_BLOCK_SIZE - RS_MAX_REDUNDANCY_BYTES)

/**
 * @brief Decode message using Reed-Solomon FEC (in-place)
 * @param *data Input/output buffer of N=K+T bytes, where K = data size and T = number of redundancy bytes
 * @param size Block size = N
 * @param *generator Generator polynomial (buffer size = T + 1)
 * @param redundancyBytes Number of redundancy bytes = T
 * @param *fixed Output number of fixed bytes
 * @return True on success, false on failure
 */
bool RsDecode(uint8_t *data, uint8_t size, uint8_t *generator, uint8_t redundancyBytes, uint8_t *fixed);

/**
 * @brief Encode message using Reed-Solomon FEC
 * @param *data Input/output buffer of N=K+T bytes, where K = data size and T = number of redundancy bytes
 * @param size Block size = N
 * @param *generator Generator polynomial (buffer size = T + 1)
 * @param redundancyBytes Number of redundancy bytes = T
 */
void RsEncode(uint8_t *data, uint8_t size, uint8_t *generator, uint8_t redundancyBytes);

/**
 * @brief Create generator polynomial for Reed-Solomon FEC
 * @param t Number of redundancy bytes = T = N - K (block size - data size)
 * @param *out Output polynomial (buffer size = T + 1), must be preallocated by caller
 * @attention Polynomial length is always 1 byte bigger than T
*/
void RsGeneratePolynomial(uint8_t t, uint8_t *out);
