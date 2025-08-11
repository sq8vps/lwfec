/*
Copyright 2023-2025 Piotr Wilkoń
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

#include "rs.h"
#include "gf.h"
#include <string.h>

//#define RS_USE_ALTERNATIVE_BM //use alternative Berlekamp-Massey implementation. Seems to be a bit faster
#define RS_USE_HORNER //use standard polynomial evalution method (Horner scheme) instead of Chien search. A bit slower

/**
 * @brief Calculates message syndromes
 * @param *rs RS instance
 * @param data Input block (length = N)
 * @param size Block size = N
 * @param out Output syndromes (length = T)
 */
static void syndromes(struct LwFecRS *rs, uint8_t *data, uint8_t size, uint8_t *out)
{
    for(uint8_t i = 0; i < rs->T; i++)
    {
        out[i] = GfPolyEval(data, size, GfPow2(i + rs->fcr));
    }
}

/**
 * @brief Calculate the error evaulator (list of erroneous positions)
 * @param *locator Error locator polynomial
 * @param locatorSize Error locator polynomial length <= T
 * @param out List of erroneous positions (error evaulator) (length = locatorSize - 1)
 * @return True on success, else the "out" buffer must be invalidated and the block is uncorrectable
 */
static bool errorEvaluator(uint8_t *locator, uint8_t locatorSize, uint8_t *out)
{
    /*
     * This function basically looks for error locator polynomial roots.
     * First implementation uses brute-force polynomial evaluation with GfPolyEval, which uses Horner's method underneath
     * Seconds implementation uses Chien search
     */
    uint8_t pos = 0;
    for(uint8_t i = 0; i < RS_BLOCK_SIZE; i++)
    {   
        #ifdef RS_USE_HORNER
        //standard evalution with Horner's method
        //evaluate at 2^i. GfPow() can be used, but taking values from exp table directly is faster
        if(GfPolyEval(locator, locatorSize, GfPow2(i % 255)) == 0) //if 2^i is the root of the locator polynomial, it determines the error position
        {
            if(pos < (locatorSize - 1))
                out[pos] = RS_BLOCK_SIZE - i - 1; //calculate error position
            else
                break;
            pos++;
        }
        #else
        //evalution with Chien search
        uint8_t lambda = 0;
        for(uint8_t j = 0; j < locatorSize; j++)
        {
            lambda ^= GfPow2((GfLog[locator[locatorSize - j - 1]] + i * j) % 255);
        }
        if(lambda == 0) 
            out[pos++] = RS_BLOCK_SIZE - i - 1; 
        #endif
    }


    if(pos != (locatorSize - 1))
        return false;

    return true;
}


/**
 * @brief Calculates the error locator polynomial
 * @param *rs RS instance
 * @param *syndromes Syndrome polynomial (length = T)
 * @param *out Output error locator buffer
 * @param outSize Error locator polynomial buffer length <= T
 * @return True if success, else the "out" buffer must be invalidated and the block is uncorrectable
 */
static bool errorLocator(struct LwFecRS *rs, uint8_t *syndromes, uint8_t *out, uint8_t *outSize)
{
    /*
     * The error locator polynomial is calculated using Berlekamp-Massey algorithm.
     * Two implementations are written here:
     * 1. directly adapted from Wikipedia (https://en.wikipedia.org/wiki/Berlekamp%E2%80%93Massey_algorithm)
     * 2. adapted/ported from Python from Wikiversity "Reed-Solomon codes for coders" 
     * (https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders#Error_correction)
     * Both implementations work just fine.
    */
#ifndef RS_USE_ALTERNATIVE_BM
    uint8_t L = 0; //number of assumed errors
    uint8_t m = 1; //number of iterations since L, B and b were updated
    uint8_t b = 1; //last discrepancy delta

    //initialize B and C polynomials with the constant term
    memset(rs->locatorBuffer.B, 0, rs->T + 1);
    rs->locatorBuffer.B[0] = 1;
    memset(rs->locatorBuffer.C, 0, rs->T + 1);
    rs->locatorBuffer.C[0] = 1;
    memset(rs->locatorBuffer.T, 0, rs->T + 1);
    memset(rs->locatorBuffer.T2, 0, rs->T + 1);

    for(uint8_t i = 0; i < rs->T; i++)
    {
        uint8_t d = syndromes[i];
        for(uint8_t j = 1; j <= L; j++)
        {
            d = GfAdd(d, GfMul(rs->locatorBuffer.C[j], syndromes[i - j])); //calculate discrepancy delta
        }
        if(d == 0) //no error
        {
            m++;
        }
        else if((L << 1) <= i)
        {
            memcpy(rs->locatorBuffer.T, rs->locatorBuffer.C, rs->T); //store C polynomial in T
            //in general, C(x)=C(x)-(d/b)*B(x)*x^m
            //first B(x)=B(x)*x^m
            //here we store polynomials as the lowest degree term first
            //multipyling it by x^m shifts the polynomial coefficient by m places right
            //so swap places starting from the last element and fill first m places with zeros
            for(uint8_t j = 0; j < (rs->T - m); j++)
            {
                rs->locatorBuffer.B[rs->T - j - 1] = rs->locatorBuffer.B[rs->T - j - m - 1];
            }
            for(uint8_t j = 0; j < m; j++)
                rs->locatorBuffer.B[j] = 0;

            //then B(x)*d/b
            GfPolyScale(rs->locatorBuffer.B, rs->T, GfMul(d, GfInv(b)), rs->locatorBuffer.B);

            //then C(x)=C(x)-B(x) (subtraction is the same as addition in GF)
            GfPolyAdd(rs->locatorBuffer.T, rs->T, rs->locatorBuffer.B, rs->T, rs->locatorBuffer.C);

            //store T polynomial in B (previous C to B)
            memcpy(rs->locatorBuffer.B, rs->locatorBuffer.T, rs->T);

            L = i + 1 - L;
            b = d;
            m = 1;
        }
        else
        {
            //the same as above, but save B in T2 first
            memcpy(rs->locatorBuffer.T, rs->locatorBuffer.C, rs->T);
            memcpy(rs->locatorBuffer.T2, rs->locatorBuffer.B, rs->T);
            for(uint8_t j = 0; j < (rs->T - m); j++)
            {
                rs->locatorBuffer.B[rs->T - j - 1] = rs->locatorBuffer.B[rs->T - j - m - 1];
            }
            for(uint8_t j = 0; j < m; j++)
                rs->locatorBuffer.B[j] = 0;
            GfPolyScale(rs->locatorBuffer.B, rs->T, GfMul(d, GfInv(b)), rs->locatorBuffer.B);
            GfPolyAdd(rs->locatorBuffer.T, rs->T, rs->locatorBuffer.B, rs->T, rs->locatorBuffer.C);
            //restore T2 to B
            memcpy(rs->locatorBuffer.B, rs->locatorBuffer.T2, rs->T);
            m++;
        }
    }

    *outSize = L + 1;
    memcpy(out, rs->locatorBuffer.C, rs->T);

    if((L * 2) > rs->T)
        return false;

    return true;
#else
    memset(rs->locatorBuffer.errLoc, 0, rs->T + 1);
    memset(rs->locatorBuffer.newLoc, 0, rs->T + 1);
    memset(rs->locatorBuffer.oldLoc, 0, rs->T + 1);
    memset(rs->locatorBuffer.tmpLoc, 0, rs->T + 1);

    uint8_t newLocLen = 0;
    uint8_t errLocLen = 1;
    uint8_t oldLocLen = 1;

    rs->locatorBuffer.errLoc[0] = 1;
    rs->locatorBuffer.oldLoc[0] = 1;

    for(uint8_t i = 0; i < rs->T; i++)
    {
        uint8_t delta = syndromes[i];
        for(uint8_t j = 1; j < errLocLen; j++)
        {
            delta = GfSub(delta, GfMul(rs->locatorBuffer.errLoc[j], syndromes[i - j]));
        }

        for(uint8_t j = 0; j < oldLocLen; j++)
        {
            rs->locatorBuffer.oldLoc[oldLocLen - j] = rs->locatorBuffer.oldLoc[oldLocLen - j - 1];
        }
        rs->locatorBuffer.oldLoc[0] = 0;

        oldLocLen++;

        if(delta != 0)
        {
            if(oldLocLen > errLocLen)
            {
                GfPolyScale(rs->locatorBuffer.oldLoc, oldLocLen, delta, rs->locatorBuffer.newLoc);
                newLocLen = oldLocLen;
                GfPolyScale(rs->locatorBuffer.errLoc, errLocLen, GfInv(delta), rs->locatorBuffer.oldLoc);
                oldLocLen = errLocLen;
                memcpy(rs->locatorBuffer.errLoc, rs->locatorBuffer.newLoc, newLocLen);
                errLocLen = newLocLen;
            }
            GfPolyScale(rs->locatorBuffer.oldLoc, oldLocLen, delta, rs->locatorBuffer.newLoc);
            newLocLen = oldLocLen;
            memcpy(rs->locatorBuffer.tmpLoc, rs->locatorBuffer.errLoc, errLocLen);
            errLocLen = GfPolyAdd(rs->locatorBuffer.tmpLoc, errLocLen, rs->locatorBuffer.newLoc, newLocLen, rs->locatorBuffer.errLoc);
        }
    }

    uint8_t index = 0;
    for(uint8_t i = 0; i < errLocLen; i++)
    {
        if((index == 0) && (rs->locatorBuffer.errLoc[i] == 0)) //drop leading zeros
            continue;

        out[index++] = rs->locatorBuffer.errLoc[i];
    }

    if(((errLocLen - 1) << 1) > rs->T)
        return false;

    *outSize = errLocLen;
    return true;
#endif
}

/**
 * @brief Calculates error magnitude (errata) polynomial and fix data
 * @param *rs RS instance
 * @param *data Input data block
 * @param size Block size = N
 * @param *syn Syndrome polynomial
 * @param *evaluator Error evaluator polynomial
 * @param errCount Number of errors (error evaulator size)
 * @return True on success, false on failure
 */
static bool fix(struct LwFecRS *rs, uint8_t *data, uint8_t size, uint8_t *syn, uint8_t *evaluator, uint8_t errCount)
{
    /*
     * This is based on Forney's algorithm.
     */

    memset(rs->fixerBuffer.locator, 0, rs->T + 1);
    memset(rs->fixerBuffer.errataEvaluator, 0, 2 * rs->T + 2);
    
    rs->fixerBuffer.locator[0] = 1; //initialize error locator to constant
    uint8_t locatorSize = 1;

    //use "errataEvaluator" as temporary variable
    for(uint8_t i = 0; i < errCount; i++)
    {
        memcpy(rs->fixerBuffer.errataEvaluator, rs->fixerBuffer.locator, rs->T);
        GfPolyMul(rs->fixerBuffer.errataEvaluator, locatorSize, (uint8_t[2]){GfPow(2, size - 1 - evaluator[i]), 1}, 2, rs->fixerBuffer.locator);
        locatorSize++;
    }

    memset(rs->fixerBuffer.errataEvaluator, 0, 2 * rs->T + 2);

    GfPolyInv(syn, rs->T);
    GfPolyMul(syn, rs->T, rs->fixerBuffer.locator, locatorSize, rs->fixerBuffer.errataEvaluator);
    for(uint8_t i = 0; i < locatorSize; i++)
    {
        rs->fixerBuffer.errataEvaluator[i] = rs->fixerBuffer.errataEvaluator[rs->T + i];
    }

    for(uint8_t i = 0; i < errCount; i++)
    {
        rs->fixerBuffer.errataPosition[i] = GfPow(2, size - 1 - evaluator[i]);
    }

    uint8_t *errLocPrimePoly = syn; //reuse
    uint8_t errLocPrimePolyLen = 0;
    for(uint8_t i = 0; i < errCount; i++)
    {
        errLocPrimePolyLen = 0;
        uint8_t errataInv = GfInv(rs->fixerBuffer.errataPosition[i]);
        for(uint8_t j = 0; j < errCount; j++)
        {
            if(j != i)
            {
                errLocPrimePoly[errLocPrimePolyLen++] = GfSub(1, GfMul(errataInv, rs->fixerBuffer.errataPosition[j]));
            }
        }
        uint8_t errLocPrime = 1;
        for(uint8_t j = 0; j < errLocPrimePolyLen; j++)
        {
            errLocPrime = GfMul(errLocPrime, errLocPrimePoly[j]);
        }

        uint8_t y = GfPolyEval(rs->fixerBuffer.errataEvaluator, locatorSize, errataInv);
        //in general y*=errataInv**(fcr-1)
        //for fcr=0 y*=errataInv**-1=errataPosition
        //for fcr=1 y does not change
        if(rs->fcr == 0)
            y = GfMul(y, rs->fixerBuffer.errataPosition[i]);
        else if(rs->fcr > 0)
            y = GfMul(y, GfPow(errataInv, rs->fcr - 1));

        if(errLocPrime == 0)
        {
            return false;
        }
        data[evaluator[i]] = GfSub(data[evaluator[i]], GfDiv(y, errLocPrime));
    }

    return true;
}


/**
 * @brief Check if syndromes are all zero, that is if the message is correct
 * @param *syndromes Syndrome polynomial
 * @param size Syndrome polynomial size (buffer length)
 * @return True if all zero
*/
static bool checkSyndromes(uint8_t *syndromes, uint8_t size)
{
    bool err = false;

    for(uint8_t i = 0; i < size; i++) //check if all syndromes are 0, if so, the message is correct
    {
        if(syndromes[i] != 0)
        {
            err = true;
        }
    }
    return !err;
}

bool RsDecode(struct LwFecRS *rs, uint8_t *data, uint8_t size, uint8_t *fixed)
{
    if((size > (RS_BLOCK_SIZE - rs->T)) || (rs->T > RS_MAX_REDUNDANCY_BYTES))
        return false;
    
    memset(rs->decodeBuffer.syndromes, 0, rs->T + 1);
    memset(rs->decodeBuffer.locator, 0, rs->T + 1);
    memset(rs->decodeBuffer.evaluator, 0, rs->T + 1);

    memmove(&data[RS_BLOCK_SIZE - rs->T], &data[size], rs->T);
    memset(&data[size], 0, RS_BLOCK_SIZE - size - rs->T);
    
    syndromes(rs, data, RS_BLOCK_SIZE, rs->decodeBuffer.syndromes); //calculate syndromes

    if(checkSyndromes(rs->decodeBuffer.syndromes, rs->T))
        return true;

    uint8_t locatorSize = 0;

    if(!errorLocator(rs, rs->decodeBuffer.syndromes, rs->decodeBuffer.locator, &locatorSize)) //calculate error locator polynomial
        return false;

    if(!errorEvaluator(rs->decodeBuffer.locator, locatorSize, rs->decodeBuffer.evaluator)) //calculate error evaulator (list of erroneous positions)
        return false;

    if(!fix(rs, data, RS_BLOCK_SIZE, 
            rs->decodeBuffer.syndromes, rs->decodeBuffer.evaluator, locatorSize - 1)) //calculate error magnitude (errata) polynomial and try to fix
        return false;

    syndromes(rs, data, RS_BLOCK_SIZE, rs->decodeBuffer.syndromes); //calculate syndromes again to check if the message has been corrected successfully

    if(checkSyndromes(rs->decodeBuffer.syndromes, rs->T))
    {
        if(NULL != fixed)
            *fixed = locatorSize - 1;
        return true;
    }
    else
        return false;
}

void RsEncode(struct LwFecRS *rs, uint8_t *data, uint8_t size)
{
    if((size > (RS_BLOCK_SIZE - rs->T)) || (rs->T > RS_MAX_REDUNDANCY_BYTES))
        return;

    memset(&data[size], 0, RS_BLOCK_SIZE - size);
    memset(rs->encodeBuffer, 0, sizeof(rs->encodeBuffer));
    memcpy(&data[size], GfPolyDiv(data, RS_BLOCK_SIZE, rs->generator, rs->T + 1, rs->encodeBuffer), rs->T);
}

void RsInit(struct LwFecRS *rs, uint8_t T, uint8_t fcr)
{
    if(T > RS_MAX_REDUNDANCY_BYTES)
        return;
    
    memset(rs->generator, 0, T + 1);
    rs->generator[0] = 1;
    for(uint8_t i = 0; i < T; i++)
    {
        memcpy(rs->initBuffer, rs->generator, i + 1);
        GfPolyMul(rs->initBuffer, i + 1, (uint8_t[2]){1, GfPow2(i + fcr)}, 2, rs->generator);
    }
    rs->T = T;
    rs->fcr = fcr;
}

#if RS_BLOCK_SIZE > 256
#error Architectural limit of RS FEC block size is 256 bytes
#endif

#if RS_MAX_DATA_SIZE <= 0
#error RS FEC parity byte count must be less than the block size
#endif