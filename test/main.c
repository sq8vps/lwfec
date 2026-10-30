
#ifdef LWFEC_TEST

#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <string.h>

#include "rs.h"

struct LwFecRS rs;

int main(int, char**)
{
    RsInit(&rs, 10, 0);
    uint8_t data[RS_BLOCK_SIZE] = "LwFEC is a lightweight FEC library for embedded systems";
    size_t size = strlen((char*)data);

    RsEncode(&rs, data, size);

    data[0] = '?';
    data[1] = '!';
    data[18] = 'j';
    data[27] = 'x';
    data[28] = 'u';

    printf("Broken data: %.*s\n", (int)size, data);

    uint8_t fixed = 0;
    if(RsDecode(&rs, data, size, &fixed))
    {
        printf("Fixed %d bytes: %.*s\n", (int)fixed, (int)size, data);
    }
    else
    {
        printf("Failed to fix data\n");
    }

    return 0;
}

#endif