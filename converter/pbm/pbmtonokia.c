/* pbmtonokia.c - convert a PBM image to Nokia Smart Messaging
   Formats (NOL, NGG, HEX)

   Copyright information is at end of file.
*/

#define _BSD_SOURCE    /* Make sure strcasecmp() is in string.h */
#include <string.h>

#include "pm_c_util.h"
#include "nstring.h"
#include "mallocvar.h"
#include "shhopt.h"
#include "pbm.h"

enum outputFormat {
    FMT_HEX_NOL,
    FMT_HEX_NGG,
    FMT_HEX_NPM,
    FMT_NOL,
    FMT_NGG
};


struct cmdlineInfo {
    /* All the information the user supplied in the command line,
       in a form easy for the program to use.
    */
    const char * inputFileName;  /* Filename of input files */
    int outputFormat;
    const char * networkCode;
    const char * txt;  /* NULL means unspecified */
};



static const char *
uppercase(const char * const subject) {

    char * buffer;

    buffer = malloc(strlen(subject) + 1);

    if (buffer == NULL)
        pm_error("Out of memory allocating buffer for uppercasing a "
                 "%u-character string", strlen(subject));
    else {
        unsigned int i;

        i = 0;
        while (subject[i]) {
            buffer[i] = TOUPPER(subject[i]);
            ++i;
        }
        buffer[i] = '\0';
    }
    return buffer;
}



static bool
ishexstring(const char * const subject) {

    bool retval;
    unsigned int i;

    retval = TRUE;  /* initial assumption */

    for (i = 0; i < strlen(subject); ++i)
        if (!ISXDIGIT(subject[i]))
            retval = FALSE;

    return retval;
}



static void
parseCommandLine(int argc, char ** argv,
                 struct cmdlineInfo * const cmdlineP) {
/*----------------------------------------------------------------------------
   Note that the file spec array we return is stored in the storage that
   was passed to us as the argv array.
-----------------------------------------------------------------------------*/
    optEntry * option_def;
        /* Instructions to optParseOptions3 on how to parse our options.
         */
    optStruct3 opt;

    unsigned int option_def_index;
    unsigned int fmtSpec, netSpec, txtSpec;
    const char * fmtOpt;
    const char * netOpt;

    MALLOCARRAY_NOFAIL(option_def, 100);

    option_def_index = 0;   /* incremented by OPTENT3 */
    OPTENT3(0, "fmt",     OPT_STRING, &fmtOpt, 
            &fmtSpec, 0);
    OPTENT3(0, "net",     OPT_STRING, &netOpt,
            &netSpec, 0);
    OPTENT3(0, "txt",     OPT_STRING, &cmdlineP->txt,
            &txtSpec, 0);

    opt.opt_table = option_def;
    opt.short_allowed = FALSE;  /* We have no short (old-fashioned) options */
    opt.allowNegNum = FALSE;  /* We have no parms that are negative numbers */

    optParseOptions3(&argc, argv, opt, sizeof(opt), 0);
        /* Uses and sets argc, argv, and some of *cmdlineP and others. */

    if (fmtSpec) {
        if (STRCASEEQ(fmtOpt, "HEX_NOL"))
            cmdlineP->outputFormat = FMT_HEX_NOL;
        else if (STRCASEEQ(fmtOpt, "HEX_NGG"))
            cmdlineP->outputFormat = FMT_HEX_NGG;
        else if (STRCASEEQ(fmtOpt, "HEX_NPM"))
            cmdlineP->outputFormat = FMT_HEX_NPM;
        else if (STRCASEEQ(fmtOpt, "NOL"))
            cmdlineP->outputFormat = FMT_NOL;
        else if (STRCASEEQ(fmtOpt, "NGG"))
            cmdlineP->outputFormat = FMT_NGG;
        else
            pm_error("-fmt option must be HEX_NOL, HEX_NGG, HEX_NPM, "
                     "NOL, or NGG.  You specified '%s'", fmtOpt);
    } else
        cmdlineP->outputFormat = FMT_HEX_NOL;

    if (netSpec) {
        if (strlen(netOpt) != 6)
            pm_error("-net option must be 6 hex digits long.  "
                     "You specified %u characters", strlen(netOpt));
        else if (!ishexstring(netOpt))
            pm_error("-net option must be hexadecimal.  You specified '%s'",
                     netOpt);
        else
            cmdlineP->networkCode = uppercase(netOpt);
    } else
        cmdlineP->networkCode = strdup("62F210");  /* German D1 net */

    if (!txtSpec)
        cmdlineP->txt = NULL;

    if (argc-1 == 0) 
        cmdlineP->inputFileName = "-";
    else if (argc-1 != 1)
        pm_error("Program takes zero or one argument (filename).  You "
                 "specified %d", argc-1);
    else
        cmdlineP->inputFileName = argv[1];
}



static void
freeCmdline(struct cmdlineInfo const cmdline) {

    strfree(cmdline.networkCode);
}



static void
convertToHexNol(bit **       const image,
                unsigned int const cols,
                unsigned int const rows,
                const char * const networkCode,
                FILE *       const ofP) {

    unsigned int row;

    /* header */
    fprintf(ofP, "06050415820000%s00%02X%02X01", networkCode, cols, rows);
    
    /* image */
    for (row = 0; row < rows; ++row) {
        unsigned int col;
        unsigned int p;
        unsigned int c;

        c = 0;

        for (p = 0, col = 0; col < cols; ++col) {
            if (image[row][col] == PBM_BLACK)
                c |= 0x80 >> p;
            if (++p == 8) {
                fprintf(ofP, "%02X",c);
                p = c = 0;
            }
        }
        if (p > 0)
            fprintf(ofP, "%02X", c);
    }
}



static void
convertToHexNgg(bit **       const image,
                unsigned int const cols,
                unsigned int const rows,
                FILE *       const ofP) {

    unsigned int row;

    /* header */
    fprintf(ofP, "0605041583000000%02X%02X01", cols, rows);

    /* image */
    for (row = 0; row < rows; ++row) {
        unsigned int col;
        unsigned int p;
        unsigned int c;

        for (p = 0, c = 0, col = 0; col < cols; ++col) {
            if (image[row][col] == PBM_BLACK)
                c |= 0x80 >> p;
            if (++p == 8) {
                fprintf(ofP, "%02X", c);
                p = c = 0;
            }
        }
        if (p > 0)
            fprintf(ofP, "%02X", c);
    }
}




static void
convertToHexNpm(bit **       const image,
                unsigned int const cols,
                unsigned int const rows,
                const char * const text,
                FILE *       const ofP) {

    unsigned int row;
    
    /* header */
    fprintf(ofP, "060504158A0000");

    /* text */
    if (text) {
        size_t const len = strlen(text);

        unsigned int it;

        fprintf(ofP, "00%04X", len);

        for (it = 0; it < len; ++it)
            fprintf(ofP, "%02X", text[it]);
    }

    /* image */
    fprintf(ofP, "02%04X00%02X%02X01", (cols * rows) / 8 + 4, cols, rows);

    for (row = 0; row < rows; ++row) {
        unsigned int col;
        unsigned int p;
        unsigned int c;

        for (p = 0, c = 0, col = 0; col < cols; ++col) {
            if (image[row][col] == PBM_BLACK)
                c |= 0x80 >> p;
            if (++p == 8) {
                fprintf(ofP, "%02X", c);
                p = c = 0;
            }
        }
        if (p > 0)
            fprintf(ofP, "%02X", c);
    }
}



static void
convertToNol(bit **       const image,
             unsigned int const cols,
             unsigned int const rows,
             FILE *       const ofP) {

    unsigned int row;
    char header[32];
    
    /* header - this is a hack */

    header[ 0] = 'N';
    header[ 1] = 'O';
    header[ 2] = 'L';
    header[ 3] = 0;
    header[ 4] = 1;
    header[ 5] = 0;
    header[ 6] = 4;
    header[ 7] = 1;
    header[ 8] = 1;
    header[ 9] = 0;
    header[10] = cols;
    header[11] = 0;
    header[12] = rows;
    header[13] = 0;
    header[14] = 1;
    header[15] = 0;
    header[16] = 1;
    header[17] = 0;
    header[18] = 0x53;
    header[19] = 0;

    fwrite(header, 20, 1, ofP);
    
    /* image */
    for (row = 0; row < rows; ++row) {
        unsigned int col;
        unsigned int p;
        unsigned int c;

        for (p = 0, c = 0, col = 0; col < cols; ++col) {
            char const output = image[row][col] == PBM_BLACK ? '1' : '0';

            putc(output, ofP);
        }
    }
}




static void
convertToNgg(bit **       const image,
             unsigned int const cols,
             unsigned int const rows,
             FILE *       const ofP) {

    unsigned int row;
    char    header[32];

    /* header - this is a hack */

    header[ 0] = 'N';
    header[ 1] = 'G';
    header[ 2] = 'G';
    header[ 3] = 0;
    header[ 4] = 1;
    header[ 5] = 0;
    header[ 6] = cols;
    header[ 7] = 0;
    header[ 8] = rows;
    header[ 9] = 0;
    header[10] = 1;
    header[11] = 0;
    header[12] = 1;
    header[13] = 0;
    header[14] = 0x4a;
    header[15] = 0;

    fwrite(header, 16, 1, ofP);
    
    /* image */

    for (row = 0; row < rows; ++row) {
        unsigned int col;
        unsigned int p;
        unsigned int c;

        for (p = 0, c = 0, col = 0; col < cols; ++col) {
            char const output = image[row][col] == PBM_BLACK ? '1' : '0';

            putc(output, ofP);
        }
    }
}



int 
main(int    argc,
     char * argv[]) {

    struct cmdlineInfo cmdline;
    FILE  * ifP;
    bit ** bits;
    int rows, cols;

    pbm_init(&argc, argv);

    parseCommandLine(argc, argv, &cmdline);

    ifP = pm_openr(cmdline.inputFileName);
    bits = pbm_readpbm(ifP, &cols, &rows);
    pm_close(ifP);

    switch (cmdline.outputFormat) {
    case FMT_HEX_NOL:
        convertToHexNol(bits, cols, rows, cmdline.networkCode, stdout);
        break;
    case FMT_HEX_NGG:
        convertToHexNgg(bits, cols, rows, stdout);
        break;
    case FMT_HEX_NPM:
        convertToHexNpm(bits, cols, rows, cmdline.txt, stdout);
        break;
    case FMT_NOL:
        convertToNol(bits, cols, rows, stdout);
        break;
    case FMT_NGG:
        convertToNgg(bits, cols, rows, stdout);
        break;
    }

freeCmdline(cmdline);

    return 0;
}



/* Copyright (C)2001 OMS Open Media System GmbH, Tim R�hsen
** <tim.ruehsen@openmediasystem.de>.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty.

  Created 2001.06.07

Notes:
  - limited to rows <= 255 and columns <= 255
  - limited to b/w graphics, not animated

Testing:
  Testing was done with SwissCom SMSC (Switzerland) and IC3S SMSC (Germany).
  The data was send with EMI/UCP protocol over TCP/IP.

  - 7.6.2001: tested with Nokia 3210: 72x14 Operator Logo
  - 7.6.2001: tested with Nokia 6210: 72x14 Operator Logo and 
              72x14 Group Graphic
*/
