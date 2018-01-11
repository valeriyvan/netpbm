/*=============================================================================
                              pamtogif
===============================================================================
  Convert a Netpbm image to GIF

  History and copyright information is at the end of the file.
=============================================================================*/

#include <assert.h>
#include <string.h>
#include <stdbool.h>

#include "pm_c_util.h"
#include "mallocvar.h"
#include "shhopt.h"
#include "nstring.h"
#include "pam.h"
#include "pammap.h"

#define MAXCMAPSIZE 256

static unsigned int const gifMaxval = 255;

static bool verbose;


typedef unsigned int StringCode;
    /* A code to be place in the GIF raster.  It represents
       a string of one or more pixels.  You interpret this in the context
       of a current code size.  The lower half of the values representable
       in the current code size represent singleton strings and the value
       is simply the value of the one pixel in the string.  The first two
       values in the upper half of the range are the clear code and EOF
       code, respectively.  The rest of the values represent multi-pixel
       strings.  The mapping between value and the sequence of pixels
       changes throughout the image.

       Ergo, this data structure must be at least BITS bits wide.
    */


struct Cmap {
    /* This is the information for the GIF colormap (aka palette). */

    struct pam pam;
        /* Gives depth and maxval for colors in color[] */
    tuple color[MAXCMAPSIZE];
        /* Maps a color index, as is found in the raster part of the
           GIF, to color.
        */
    unsigned int cmapSize;
        /* Number of entries in the GIF colormap.  I.e. number of colors
           in the image, plus possibly one fake transparency color.
        */
    bool haveTransparent;
        /* The colormap contains an entry for transparent pixels */
    unsigned int transparent;
        /* color index number in GIF palette of the color that is to
           be transparent.

           Meaningful only if 'haveTransparent' is true.
        */
    tuplehash tuplehash;
        /* A hash table to translate color to GIF colormap index. */
};

struct CmdlineInfo {
    /* All the information the user supplied in the command line,
       in a form easy for the program to use.
    */
    const char *input_filespec; /* Filespec of input file */
    const char *alphacolor;     /* -alphacolor option value or default */
    unsigned int interlace;     /* -interlace option value */
    unsigned int sort;          /* -sort option value */
    const char *mapfile;        /* -mapfile option value.  NULL if none. */
    const char *transparent;    /* -transparent option value.  NULL if none. */
    const char *comment;        /* -comment option value; NULL if none */
    unsigned int nolzw;         /* -nolzw option */
    float aspect;               /* -aspect option value (the ratio).  */
    unsigned int verbose;
};





static unsigned int
pamAlphaPlane(struct pam * const pamP) {

    unsigned int alphaPlane;

    if (streq(pamP->tuple_type, "RGB_ALPHA"))
        alphaPlane = 3;
    else if (streq(pamP->tuple_type, "GRAYSCALE_ALPHA"))
        alphaPlane = 2;
    else if (streq(pamP->tuple_type, "BLACKANDWHITE_ALPHA"))
        alphaPlane = 2;
    else
        alphaPlane = 0;

    if (alphaPlane >= pamP->depth)
        pm_error("Tuple type is '%s', but depth (%u) is less than %u",
                 pamP->tuple_type, pamP->depth, alphaPlane + 1);

    return alphaPlane;
}



static void
parseCommandLine(int argc, char ** argv,
                 struct CmdlineInfo * const cmdlineP) {
/*----------------------------------------------------------------------------
   Parse the program arguments (given by argc and argv) into a form
   the program can deal with more easily -- a cmdline_info structure.
   If the syntax is invalid, issue a message and exit the program via
   pm_error().

   Note that the file spec array we return is stored in the storage that
   was passed to us as the argv array.
-----------------------------------------------------------------------------*/
    optEntry * option_def;  /* malloc'ed */
    optStruct3 opt;  /* set by OPTENT3 */
    unsigned int option_def_index;

    unsigned int aspectSpec;

    MALLOCARRAY_NOFAIL(option_def, 100);

    option_def_index = 0;   /* incremented by OPTENT3 */
    OPTENT3(0,   "interlace",   OPT_FLAG,
            NULL,                       &cmdlineP->interlace, 0);
    OPTENT3(0,   "sort",        OPT_FLAG,
            NULL,                       &cmdlineP->sort, 0);
    OPTENT3(0,   "nolzw",       OPT_FLAG,
            NULL,                       &cmdlineP->nolzw, 0);
    OPTENT3(0,   "mapfile",     OPT_STRING,
            &cmdlineP->mapfile,        NULL, 0);
    OPTENT3(0,   "transparent", OPT_STRING,
            &cmdlineP->transparent,    NULL, 0);
    OPTENT3(0,   "comment",     OPT_STRING,
            &cmdlineP->comment,        NULL, 0);
    OPTENT3(0,   "alphacolor",  OPT_STRING,
            &cmdlineP->alphacolor,     NULL, 0);
    OPTENT3(0,   "aspect",      OPT_FLOAT,
            &cmdlineP->aspect,         &aspectSpec, 0);
    OPTENT3(0,   "verbose",     OPT_FLAG,
            NULL,                      &cmdlineP->verbose, 0);

    /* Set the defaults */
    cmdlineP->mapfile = NULL;
    cmdlineP->transparent = NULL;  /* no transparency */
    cmdlineP->comment = NULL;      /* no comment */
    cmdlineP->alphacolor = "rgb:0/0/0";
        /* We could say "black" here, but then we depend on the color names
           database existing.
        */

    opt.opt_table = option_def;
    opt.short_allowed = FALSE;  /* We have no short (old-fashioned) options */
    opt.allowNegNum = FALSE;  /* We have no parms that are negative numbers */

    pm_optParseOptions3(&argc, argv, opt, sizeof(opt), 0);
        /* Uses and sets argc, argv, and some of *cmdlineP and others. */

    if (argc-1 == 0)
        cmdlineP->input_filespec = "-";
    else if (argc-1 != 1)
        pm_error("Program takes zero or one argument (filename).  You "
                 "specified %d", argc-1);
    else
        cmdlineP->input_filespec = argv[1];

    if (aspectSpec) {
        if (cmdlineP->aspect < 0.25  || cmdlineP->aspect > 4.21875)
            pm_error("Invalid -aspect value: %f.  "
                     "GIF allows only the range 0.25-4.0 .",
                     cmdlineP->aspect);
        else if (cmdlineP->aspect > 4.0)
            pm_message("Warning: "
                       "You specified an aspect ratio over 4.0: %f.  "
                       "This will result in an invalid GIF.",
                       cmdlineP->aspect);
    } else
        cmdlineP->aspect = 1.0;
}



/*
 * Write out a word to the GIF file
 */
static void
Putword(int const w, FILE * const fp) {

    fputc( w & 0xff, fp );
    fputc( (w / 256) & 0xff, fp );
}



static int
closestColor(tuple         const color,
             struct pam *  const pamP,
             struct Cmap * const cmapP) {
/*----------------------------------------------------------------------------
   Return the colormap index of the color in the colormap *cmapP
   that is closest to the color 'color', whose format is specified by
   *pamP.

   Also add 'color' to the colormap hash, with the colormap index we
   are returning.  Caller must ensure that the color is not already in
   there.
-----------------------------------------------------------------------------*/
    unsigned int const nComp = pamP->depth >= 3 ? 3 : 1;
        /* Number of color components (not alpha) in 'color' */

    unsigned int i;
    unsigned int imin, dmin;
    int fits;

    dmin = UINT_MAX;
    imin = 0;
    for (i = 0; i < cmapP->cmapSize; ++i) {
        unsigned int distance;
        unsigned int plane;

        for (distance = 0, plane = 0; plane < nComp; ++plane)
            /* Divide by 4 is to avoid arithmetic overflow */
            distance += SQR(color[plane] - cmapP->color[i][plane]) / 4;

        if (distance < dmin) {
            dmin = distance;
            imin = i;
        }
    }
    pnm_addtotuplehash(pamP, cmapP->tuplehash, color, imin, &fits);

    return imin;
}



enum pass {MULT8PLUS0, MULT8PLUS4, MULT4PLUS2, MULT2PLUS1};


typedef struct {
    struct pam pam;
        /* Description of input file/image.  The position of the file
           is also part of the state of this rowReader.
        */
    pm_filepos rasterPos;
        /* Position in file fileP of the start of the raster */
    bool interlace;
        /* We're accessing the image in interlace fashion */
    bool eof;
        /* The image is at EOF (we have returned all of the rows) */
    unsigned int nextRow;
        /* Number of row to which input file is positioned;
           meaningless if 'eof'.
        */
    enum pass pass;
        /* The interlace pass.  Undefined if !interlace */
    tuple * discardBuffer;
        /* A bitbucket for rows we read in order to advance the file
           position.
        */
} RowReader;



static RowReader *
rowReader_create(struct pam * const pamP,
                 pm_filepos   const rasterPos,
                 bool         const interlace) {

    RowReader * rdrP;

    MALLOCVAR_NOFAIL(rdrP);

    rdrP->pam         = *pamP;
    rdrP->rasterPos   = rasterPos;
    rdrP->interlace   = interlace;
    rdrP->eof         = FALSE;
    rdrP->pass        = MULT8PLUS0;

    pm_seek2(rdrP->pam.file, &rasterPos, sizeof(rasterPos));
    rdrP->nextRow = 0;

    rdrP->discardBuffer = pnm_allocpamrow(&rdrP->pam);

    return rdrP;
}



static void
rowReader_destroy(RowReader * const rdrP) {

    pnm_freepamrow(rdrP->discardBuffer);
    free(rdrP);
}



static void
rowReaderSkipRows(RowReader *  const rdrP,
                  unsigned int const rowCount,
                  bool *       const eofP) {
/*----------------------------------------------------------------------------
   Skip over the next 'rowCount' rows of the input file.

   Iff there aren't at least 'rowCount' rows left, return *eofP == TRUE.
-----------------------------------------------------------------------------*/
    if (rdrP->nextRow + rowCount >= rdrP->pam.height)
        *eofP = TRUE;
    else {
        /* This could be made faster if need be by adding a libnetpbm
           row skip function.  Except with the plain formats, that could
           just compute the next row position and fseek() to it.
           pnm_readpamrow() with NULL for the output pointer would be a
           good interface for a row skip function.
        */
        unsigned int i;

        *eofP = FALSE;

        for (i = 0; i < rowCount; ++i)
            pnm_readpamrow(&rdrP->pam, rdrP->discardBuffer);

        rdrP->nextRow += rowCount;
    }
}



static void
rowReaderGotoNextInterlaceRow(RowReader * const rdrP) {
/*----------------------------------------------------------------------------
  Position reader to the next row in the interlace pattern, assuming it
  is now positioned immediately after the current row.
-----------------------------------------------------------------------------*/
    bool endOfPass;

    /* There are 4 passes:
       MULT8PLUS0: Rows 0, 8, 16, 24, 32, etc.
       MULT8PLUS4: Rows 4, 12, 20, 28, etc.
       MULT4PLUS2: Rows 2, 6, 10, 14, etc.
       MULT2PLUS1: Rows 1, 3, 5, 7, 9, etc.
    */

    switch (rdrP->pass) {
    case MULT8PLUS0:
        rowReaderSkipRows(rdrP, 7, &endOfPass);
        break;
    case MULT8PLUS4:
        rowReaderSkipRows(rdrP, 7, &endOfPass);
        break;
    case MULT4PLUS2:
        rowReaderSkipRows(rdrP, 3, &endOfPass);
        break;
    case MULT2PLUS1:
        rowReaderSkipRows(rdrP, 1, &endOfPass);
        break;
    }

    /* Note that if there are more than 4 rows, the sequence of passes
       is sequential, but when there are fewer than 4, reading may skip
       e.g. from MULT8PLUS0 to MULT4PLUS2.
    */
    while (endOfPass && !rdrP->eof) {
        pm_seek2(rdrP->pam.file, &rdrP->rasterPos, sizeof(rdrP->rasterPos));
        rdrP->nextRow = 0;

        switch (rdrP->pass) {
        case MULT8PLUS0:
            rdrP->pass = MULT8PLUS4;
            rowReaderSkipRows(rdrP, 4, &endOfPass);
            break;
        case MULT8PLUS4:
            rdrP->pass = MULT4PLUS2;
            rowReaderSkipRows(rdrP, 2, &endOfPass);
            break;
        case MULT4PLUS2:
            rdrP->pass = MULT2PLUS1;
            rowReaderSkipRows(rdrP, 1, &endOfPass);
            break;
        case MULT2PLUS1:
            rdrP->eof = TRUE;
            break;
        }
    }
}




static void
rowReaderGotoNextStraightRow(RowReader * const rdrP) {
/*----------------------------------------------------------------------------
  Position reader to the next row in a straight, non-interlace
  pattern, assuming the file is now positioned immediately after the
  current row.

  This is trivial, since the next row _is_ immediately after the current
  row, except in the case that there are no more rows.
-----------------------------------------------------------------------------*/
    if (rdrP->nextRow >= rdrP->pam.height)
        rdrP->eof = TRUE;
}



static void
rowReader_read(RowReader * const rdrP,
               tuple *     const tuplerow) {

    if (rdrP->eof)
        pm_error("INTERNAL ERROR: rowReader attempted to read beyond end "
                 "of image");

    pnm_readpamrow(&rdrP->pam, tuplerow);
    ++rdrP->nextRow;

    if (rdrP->interlace)
        rowReaderGotoNextInterlaceRow(rdrP);
    else
        rowReaderGotoNextStraightRow(rdrP);
}



static unsigned int
gifPixel(struct pam *   const pamP,
         tuple          const tuple,
         unsigned int   const alphaPlane,
         sample         const alphaThreshold,
         struct Cmap *  const cmapP) {
/*----------------------------------------------------------------------------
   Return as *colorIndexP the colormap index of the tuple 'tuple',
   whose format is described by *pamP, using colormap *cmapP.

   'alphaThreshold' is the alpha level below which we consider a
   pixel transparent for GIF purposes.
-----------------------------------------------------------------------------*/
    int colorIndex;

    if (alphaPlane && tuple[alphaPlane] < alphaThreshold &&
        cmapP->haveTransparent)
        colorIndex = cmapP->transparent;
    else {
        int found;

        pnm_lookuptuple(pamP, cmapP->tuplehash, tuple,
                        &found, &colorIndex);

        if (!found)
            colorIndex = closestColor(tuple, pamP, cmapP);
    }
    assert(colorIndex >= 0);
    return (unsigned int) colorIndex;
}



static void
writeTransparentColorIndexExtension(FILE *       const ofP,
                                    unsigned int const transColorIndex) {
/*----------------------------------------------------------------------------
   Write out extension for transparent color index.
-----------------------------------------------------------------------------*/
    fputc('!',  ofP);
    fputc(0xf9, ofP);
    fputc(4,    ofP);
    fputc(1,    ofP);
    fputc(0,    ofP);
    fputc(0,    ofP);
    fputc(transColorIndex, ofP);
    fputc(0,    ofP);
}



static void
writeCommentExtension(FILE * const ofP,
                      char   const comment[]) {
/*----------------------------------------------------------------------------
   Write out extension for a comment
-----------------------------------------------------------------------------*/
    unsigned int const maxSegmentSize = 255;

    const char * segment;

    fputc('!',  ofP);   /* Identifies an extension */
    fputc(0xfe, ofP);   /* Identifies a comment */

    /* Write it out in segments no longer than 255 characters */
    for (segment = &comment[0];
         segment < comment + strlen(comment);
         segment += maxSegmentSize) {

        unsigned int const lengthThisSegment =
            MIN(maxSegmentSize, strlen(segment));

        fputc(lengthThisSegment, ofP);

        fwrite(segment, 1, lengthThisSegment, ofP);
    }

    fputc(0, ofP);   /* No more comment blocks in this extension */
}



/***************************************************************************
 *
 *  GIFCOMPR.C       - GIF Image compression routines
 *
 *  Lempel-Ziv compression based on 'compress'.  GIF modifications by
 *  David Rowley (mgardi@watdcsu.waterloo.edu)
 *
 ***************************************************************************/

/*
 * General DEFINEs
 */

#define BITS    12

/*
 *
 * GIF Image compression - modified 'compress'
 *
 * Based on: compress.c - File compression ala IEEE Computer, June 1984.
 *
 * By Authors:  Spencer W. Thomas       (decvax!harpo!utah-cs!utah-gr!thomas)
 *              Jim McKie               (decvax!mcvax!jim)
 *              Steve Davies            (decvax!vax135!petsd!peora!srd)
 *              Ken Turkowski           (decvax!decwrl!turtlevax!ken)
 *              James A. Woods          (decvax!ihnp4!ames!jaw)
 *              Joe Orost               (decvax!vax135!petsd!joe)
 *
 */


static StringCode const maxCodeLimitLzw = (StringCode)1 << BITS;
       /* One beyond the largest string code that can exist in GIF */
       /* Used only in assertions  */


struct HashTableEntry {
    /* This is an entry in the string table, which is a hash table.  It says
       that the string code 'combinedString' represents the string which is
       the single pixel 'additionalPixel' appended to 'baseString', where
       'baseString' may represent a multi-pixel string.
    */
    bool present;
        /* There is an entry here.  Following members are meaningless if
           not.
        */
    StringCode baseString;
    StringCode additionalPixel;
    StringCode combinedString;
};



/***************************************************************************
*                          BYTE OUTPUTTER
***************************************************************************/

typedef struct {
    FILE * fileP;  /* The file to which to output */
    unsigned int count;
        /* Number of bytes so far in the current data block */
    unsigned char buffer[256];
        /* The current data block, under construction */
} ByteBuffer;



static ByteBuffer *
byteBuffer_create(FILE * const fileP) {

    ByteBuffer * byteBufferP;

    MALLOCVAR_NOFAIL(byteBufferP);

    byteBufferP->fileP = fileP;
    byteBufferP->count = 0;

    return byteBufferP;
}



static void
byteBuffer_destroy(ByteBuffer * const byteBufferP) {

    free(byteBufferP);
}



static void
byteBuffer_flush(ByteBuffer * const byteBufferP) {
/*----------------------------------------------------------------------------
   Write the current data block to the output file, then reset the current
   data block to empty.
-----------------------------------------------------------------------------*/
    if (byteBufferP->count > 0 ) {
        if (verbose)
            pm_message("Writing %u byte block", byteBufferP->count);
        fputc(byteBufferP->count, byteBufferP->fileP);
        fwrite(byteBufferP->buffer, 1, byteBufferP->count, byteBufferP->fileP);
        byteBufferP->count = 0;
    }
}



static void
byteBuffer_flushFile(ByteBuffer * const byteBufferP) {

    fflush(byteBufferP->fileP);

    if (ferror(byteBufferP->fileP))
        pm_error("error writing output file");
}



static void
byteBuffer_out(ByteBuffer *  const byteBufferP,
               unsigned char const c) {
/*----------------------------------------------------------------------------
  Add a byte to the end of the current data block, and if it is now 255
  characters, flush the data block to the output file.
-----------------------------------------------------------------------------*/
    byteBufferP->buffer[byteBufferP->count++] = c;
    if (byteBufferP->count >= 255)
        byteBuffer_flush(byteBufferP);
}



/***************************************************************************
*                          GIF CODE OUTPUTTER
***************************************************************************/

typedef struct {
    ByteBuffer * byteBufferP;
    unsigned int initBits;
    unsigned int nBits;
        /* Number of bits to put in output for each code */
    StringCode maxCode;                  /* maximum code, given n_bits */
    StringCode maxCodeLimit;
        /* LZW: One beyond the largest string code that can exist in GIF.
           Uncompressed: a ceiling to prevent code size from ratcheting up.
           In either case, output code never reaches this value.
        */
    unsigned long curAccum;
    int curBits;
    unsigned int codeCount;
        /* Number of codes that have been output to this buffer (doesn't
           matter if they have gone out the other side yet or not) since
           the last flush (or ever, if no last flush).  The main use of this
           is debugging -- when something fails, you can see in a debugger
           where in the image it was, then set a trap for there.
        */
} CodeBuffer;



static CodeBuffer *
codeBuffer_create(FILE *       const ofP,
                  unsigned int const initBits,
                  bool         const lzw) {

    CodeBuffer * codeBufferP;

    MALLOCVAR_NOFAIL(codeBufferP);

    codeBufferP->initBits    = initBits;
    codeBufferP->nBits       = codeBufferP->initBits;
    codeBufferP->maxCode     = (1 << codeBufferP->nBits) - 1;
    codeBufferP->maxCodeLimit = lzw ?
        (StringCode)1 << BITS : (StringCode) (1 << codeBufferP->nBits) - 1;
    codeBufferP->byteBufferP = byteBuffer_create(ofP);
    codeBufferP->curAccum    = 0;
    codeBufferP->curBits     = 0;
    codeBufferP->codeCount   = 0;

    return codeBufferP;
}



static void
codeBuffer_destroy(CodeBuffer * const codeBufferP) {

    byteBuffer_destroy(codeBufferP->byteBufferP);

    free(codeBufferP);
}



static void
codeBuffer_resetCodeSize(CodeBuffer * const codeBufferP) {

    codeBufferP->nBits = codeBufferP->initBits;

    assert(codeBufferP->nBits <= BITS);

    codeBufferP->maxCode = (1 << codeBufferP->nBits) - 1;
}



static void
codeBuffer_increaseCodeSize(CodeBuffer * const codeBufferP) {

    ++codeBufferP->nBits;

    assert(codeBufferP->nBits <= BITS);

    codeBufferP->maxCode = (1 << codeBufferP->nBits) - 1;
}

static void
codeBuffer_output(CodeBuffer * const codeBufferP,
                  StringCode   const code) {
/*----------------------------------------------------------------------------
   Output one GIF code to the file, through the code buffer.

   The code is represented as N bits in the file -- the lower
   N bits of 'code'.  N is a the current code size of *codeBufferP.

   Id 'code' is the maximum possible code for the current code size
   for *codeBufferP, increase that code size (unless it's already
   maxed out).
-----------------------------------------------------------------------------*/
    assert (code <= codeBufferP->maxCode);

    codeBufferP->curAccum &= (1 << codeBufferP->curBits) - 1;

    if (codeBufferP->curBits > 0)
        codeBufferP->curAccum |= ((unsigned long)code << codeBufferP->curBits);
    else
        codeBufferP->curAccum = code;

    codeBufferP->curBits += codeBufferP->nBits;

    while (codeBufferP->curBits >= 8) {
        byteBuffer_out(codeBufferP->byteBufferP,
                       codeBufferP->curAccum & 0xff);
        codeBufferP->curAccum >>= 8;
        codeBufferP->curBits -= 8;
    }

    ++codeBufferP->codeCount;
}



static void
codeBuffer_flush(CodeBuffer * const codeBufferP) {

    /* Output the possible partial byte in the buffer */

    if (codeBufferP->curBits > 0) {
        byteBuffer_out(codeBufferP->byteBufferP,
                       codeBufferP->curAccum & 0xff);
        codeBufferP->curBits = 0;
    }
    byteBuffer_flush(codeBufferP->byteBufferP);

    byteBuffer_flushFile(codeBufferP->byteBufferP);

    if (verbose)
        pm_message("%u strings of pixels written to file",
                   codeBufferP->codeCount);
    codeBufferP->codeCount = 0;
}



typedef struct {
    CodeBuffer * codeBufferP;
        /* The place to which we write our string codes.

           Constant.
        */
    bool lzw;
        /* We're actually doing LZW compression.  False means we follow
           the algorithm enough tht an LZW decompressor will recover the
           proper data, but always using one code per pixel, and therefore
           not effecting any compression and not using the LZW patent.
        */
    unsigned int hsize;
        /* The number of slots in the hash table.  This variable to
           enhance overall performance by reducing memory use when
           encoding smaller gifs.
         */

    unsigned int hshift;
        /* This is how many bits we shift left a string code in forming the
           primary hash of the concatenation of that string with another.
           Constant.
        */

    /* Codes less than 'clearCode' are singleton pixel codes - each
       represents the pixel value equal to it.

       Codes greater than 'eofCode' are multipixel string codes.  Each
       represents a string of pixels that is defined by the preceding
       stream.
    */
    StringCode clearCode;
        /* The code in an LZW stream that means to clear the string
           dictionary and start fresh.

           Constant.
        */
    StringCode eofCode;
        /* The code in an LZW stream that means there's no more coming

           Constant.
        */
    StringCode initCodeLimit;
        /* The value of 'codeLimit' at the start of a block.

           Constant.
        */

    StringCode codeLimit;
        /* One beyond the maximum code possible with the current code
           size.
        */

    struct HashTableEntry * hashTable;
    StringCode nextUnusedCode;
        /* Numerically next code available to assign to a a multi-pixel
           string.  Note that codes for multi-pixel strings are in the
           upper half of the range of codes, always greater than
           'clearCode'.
        */

    StringCode stringSoFar;
        /* The code for the string we have built so far.  This code indicates
           one or more pixels that we have encoded but not yet output
           because we're hoping to match an even longer string.

           Valid only when 'buildingString' is true.

           In the non-lzw case the single pixel to output.
        */
    bool buildingString;
        /* We are in the middle of building a string; 'stringSoFar' describes
           the pixels in it so far.  The only time this is false is at the
           very beginning of the stream.

           Ignored in the non-lzw case.
        */
} LzwCompressor;




static unsigned int
nSignificantBits(unsigned int const arg){

#if HAVE_GCC_BITCOUNT

    return (arg == 0) ? 0 : 8 * sizeof(unsigned int) - __builtin_clz(arg);

#else

    unsigned int i = 0;
    while (arg >> i != 0)
        ++i;

    return i;
#endif
}



static LzwCompressor *
lzw_create(FILE *       const ofP,
           unsigned int const initBits,
           bool         const lzw,
           unsigned int const pixelCount) {

    unsigned int const hsizeTable[] = {257, 521, 1031, 2053, 4099, 5003};
    /* If the image has 4096 or fewer pixels we use prime numbers slightly
       above powers of two between 8 and 12.  In this case the hash table
       never fills up; clear code is never emitted.

       Above that we use a table with 4096 slots plus 20% extra.
       When this is not enough the clear code is emitted.
       Because of the extra 20% the table itself never fills up.

       lzw.hsize and lzw.hshift stay constant through the image.

       Variable hsize is a performance enhancement based on the fact that
       the encoder never needs more codes than the number of pixels in
       the image.  Typically, the ratio of pixels to codes is around
       10:1 to 20:1.

       Logic works with fixed values lzw.hsize=5003 and t=13.
    */

    LzwCompressor * lzwP;

    MALLOCVAR_NOFAIL(lzwP);

    /* Constants */
    lzwP->lzw = lzw;

    lzwP->clearCode     = 1 << (initBits - 1);
    lzwP->eofCode       = lzwP->clearCode + 1;
    lzwP->initCodeLimit = 1 << initBits;

    if (lzw) {
        unsigned int const t =
            MIN(13, MAX(8, nSignificantBits(pixelCount +lzwP->eofCode - 2)));
            /* Index into hsizeTable */

        lzwP->hsize = hsizeTable[t-8];

        lzwP->hshift = (t == 13 ? 12 : t) - nSignificantBits(MAXCMAPSIZE-1);

        MALLOCARRAY(lzwP->hashTable, lzwP->hsize);

        if (lzwP->hashTable == NULL)
            pm_error("Couldn't get memory for %u-entry hash table.",
                     lzwP->hsize);
    } else {
        /* No LZW compression.  We don't need a stringcode hash table */
        lzwP->hashTable = NULL;
        lzwP->hsize     = 0;
    }

    lzwP->buildingString = FALSE;

    lzwP->codeBufferP = codeBuffer_create(ofP, initBits, lzw);

    return lzwP;
}



static void
lzw_destroy(LzwCompressor * const lzwP) {

    codeBuffer_destroy(lzwP->codeBufferP);

    free(lzwP->hashTable);

    free(lzwP);
}



static void
lzwHashClear(LzwCompressor * const lzwP) {

    /* Empty the code table */

    unsigned int i;

    for (i = 0; i < lzwP->hsize; ++i)
        lzwP->hashTable[i].present = false;

    lzwP->nextUnusedCode = lzwP->clearCode + 2;
}



static void
lzw_clearBlock(LzwCompressor * const lzwP) {
/*----------------------------------------------------------------------------
  Insert a string table clear in the stream.  Clear our table and set it up to
  start building again, and emit the code to tell the decoder we're doing it
  so he can do the same.
-----------------------------------------------------------------------------*/
    lzwHashClear(lzwP);

    codeBuffer_output(lzwP->codeBufferP, lzwP->clearCode);

    codeBuffer_resetCodeSize(lzwP->codeBufferP);

    lzwP->codeLimit = lzwP->initCodeLimit;
}



static void
lzwAdjustCodeSize(LzwCompressor * const lzwP,
                  StringCode      const newCode) {
/*----------------------------------------------------------------------------
   Assuming we just defined code 'newCode', increase the code size as
   required so that this code fits.

   The decompressor is mimicking our assignment of that code, so knows that
   we are making this adjustment, so expects codes of the new size.
-----------------------------------------------------------------------------*/
    assert(newCode <= lzwP->codeLimit);

    if (newCode == lzwP->codeLimit) {
        lzwP->codeLimit *= 2;
        codeBuffer_increaseCodeSize(lzwP->codeBufferP);

        assert(lzwP->codeLimit <= maxCodeLimitLzw);
    }
}



static void
lzwOutputCurrentString(LzwCompressor * const lzwP) {
/*----------------------------------------------------------------------------
   Put a code for the currently built-up string in the output stream.

   Doing this causes a new string code to be defined (code is
   lzwP->nextUnusedCode), so Caller must add that to the hash.  If
   that code's size is beyond the overall limit, we reset the hash
   (which means future codes will start back at the minimum size) and
   put a clear code in the stream to tell the decompressor to do the
   same.  So Caller must add it to the hash _before_ calling us.

   Note that in the non-compressing case, the overall limit is small
   enough to prevent us from ever defining string codes; we'll always
   reset the hash.

   There's an odd case that always screws up any attempt to make this
   code cleaner: At the end of the LZW stream, you have to output the
   code for the final string even though you don't have a following
   pixel that would make a longer string.  So there's nothing to add
   to the hash table and no point in allocating a new string code.
   But the decompressor doesn't know that we're done, so he allocates
   the next string code and may therefore increase his code length.
   If we don't do the same, we will write our one last code -- the EOF
   code -- in a code length smaller than what the decompressor is
   expecting, and he will have a premature end of stream.

   So this subroutine does run for that final code flush and does some
   of the motions of defining a new string code, but this subroutine
   can't update the hash because in that particular case, there's
   nothing to add.
-----------------------------------------------------------------------------*/
    codeBuffer_output(lzwP->codeBufferP, lzwP->stringSoFar);
    if (lzwP->nextUnusedCode < lzwP->codeBufferP->maxCodeLimit) {
        /* Allocate the code for the extended string, which Caller
           should have already put in the hash so he can use it in the
           future.  Decompressor knows when it sees the code output
           above to define a string on its end too, using the same
           string code we do.
        */
        StringCode const newCode = lzwP->nextUnusedCode++;

        /* This code may be too big to fit in the current code size, in
           which case we have to increase the code size (and decompressor
           will do the same).
        */
        lzwAdjustCodeSize(lzwP, newCode);
    } else {
        /* Forget all the strings so far; start building again; tell
           decompressor to do the same.
        */
        lzw_clearBlock(lzwP);
    }
}



static void
lzw_flush(LzwCompressor * const lzwP) {

    if (lzwP->lzw)
        lzwOutputCurrentString(lzwP);
        /* Put out the code for the final string. */

    codeBuffer_output(lzwP->codeBufferP, lzwP->eofCode);

    codeBuffer_flush(lzwP->codeBufferP);
}



static unsigned int
primaryHash(StringCode   const baseString,
            StringCode   const additionalPixel,
            unsigned int const hshift) {

    unsigned int hash;

    assert(baseString < maxCodeLimitLzw);
    assert(additionalPixel < MAXCMAPSIZE);

    hash = (additionalPixel << hshift) ^ baseString;

    return hash;
}



static void
lookupInHash(LzwCompressor *  const lzwP,
             unsigned int     const gifPixel,
             bool *           const foundP,
             StringCode *     const codeP,
             unsigned int *   const hashP) {

    int disp;
        /* secondary hash stride (after G. Knott) */
    int i;
        /* Index into hash table */

    i = primaryHash(lzwP->stringSoFar, gifPixel, lzwP->hshift);
    disp = (i == 0) ? 1 : lzwP->hsize - i;

    while (lzwP->hashTable[i].present &&
           (lzwP->hashTable[i].baseString != lzwP->stringSoFar ||
            lzwP->hashTable[i].additionalPixel != gifPixel)) {
        i -= disp;
        if (i < 0)
            i += lzwP->hsize;
    }

    if (lzwP->hashTable[i].present) {
        /* Found fcode in hash table */
        *foundP = true;
        *codeP = lzwP->hashTable[i].combinedString;
    } else {
        /* Found where it _should_ be (but it's not) with primary hash */
        *foundP = false;
        *hashP = i;
    }
}



static void
lzw_encodePixel(LzwCompressor * const lzwP,
                unsigned int    const gifPixel) {

    assert(gifPixel < 256);

    if (!lzwP->buildingString) {
        /* Start a new string with just this pixel */
        lzwP->stringSoFar = gifPixel;
        lzwP->buildingString = true;
    } else {
        bool found;
            /* There's a code for the current string in the string table */
        StringCode code;
            /* Existing code for the current string in the string table, if
               any
            */
        unsigned int hash;
            /* Index into hash table where the entry for the new string code
               should go; meaningless if we don't need a new string code
               (because there's already one in the hash table for the
               current string)
            */

        lookupInHash(lzwP, gifPixel, &found, &code, &hash);

        if (found)
            /* With this new pixel, it is still a known string; 'code' is
               its code
            */
            lzwP->stringSoFar = code;
        else {
            /* Make up a new string code for the current string (for use if we
               see that string later in the stream), and add it to the string
               table.
            */

            lzwP->hashTable[hash].present         = true;
            lzwP->hashTable[hash].baseString      = lzwP->stringSoFar;
            lzwP->hashTable[hash].additionalPixel = gifPixel;
            lzwP->hashTable[hash].combinedString  = lzwP->nextUnusedCode;

            /* Output the code for the known prefix of the string, thus
               defining a new string code for possible later use.  Warning:
               lzwOutputCurrentString() does more than you think.
            */

            lzwOutputCurrentString(lzwP);

            /* This singleton pixel starts the next string */
            lzwP->stringSoFar = gifPixel;
        }
    }
}



/*
 * Algorithm:  use open addressing double hashing (no chaining) on the
 * prefix code / next character combination.  We do a variant of Knuth's
 * algorithm D (vol. 3, sec. 6.4) along with G. Knott's relatively-prime
 * secondary probe.  Here, the modular division first probe is gives way
 * to a faster exclusive-or manipulation.  Also do block compression with
 * an adaptive reset, whereby the code table is cleared when the compression
 * ratio decreases, but after the table fills.  The variable-length output
 * codes are re-sized at this point, and a special CLEAR code is generated
 * for the decompressor.  Late addition:  construct the table according to
 * file size for noticeable speed improvement on small files.  Please direct
 * questions about this implementation to ames!jaw.
 */

static void
writePixelUncompressed(LzwCompressor * const lzwP,
                       unsigned int    const gifPixel) {

    lzwP->stringSoFar = gifPixel;
    lzwOutputCurrentString(lzwP);

}

static void
writeRaster(struct pam *  const pamP,
            RowReader *   const rowReaderP,
            unsigned int  const alphaPlane,
            unsigned int  const alphaThreshold,
            struct Cmap * const cmapP,
            unsigned int  const initBits,
            FILE *        const ofP,
            bool          const lzw) {
/*----------------------------------------------------------------------------
   Write the raster to file 'ofP'.

   Get the raster to write from 'rowReaderP', which gives tuples whose
   format is described by 'pamP'.

   Use the colormap 'cmapP' to generate the raster ('rowReaderP' gives
   pixel values as RGB samples; the GIF raster is colormap indices).

   Write the raster using LZW compression, or uncompressed depending
   on 'lzw'.
-----------------------------------------------------------------------------*/
    LzwCompressor * lzwP;
    tuple * tuplerow;
    unsigned int nRowsDone;
        /* Number of rows we have read so far from the the input (the
           last of which is the one we're working on now).  Note that
           in case of interlace, this is not the same thing as the row
           number of the current row.
        */

    lzwP = lzw_create(ofP, initBits, lzw, pamP->height * pamP->width);

    tuplerow = pnm_allocpamrow(pamP);

    lzw_clearBlock(lzwP);

    nRowsDone = 0;

    while (nRowsDone < pamP->height) {
        unsigned int col;

        rowReader_read(rowReaderP, tuplerow);

        for (col = 0; col < pamP->width; ++col) {
            unsigned int const colorIndex =
                gifPixel(pamP, tuplerow[col], alphaPlane, alphaThreshold,
                         cmapP);

                /* The value for the pixel in the GIF image.  I.e. the colormap
                   index.
                */
            if (lzw)
                lzw_encodePixel(lzwP, colorIndex);
            else
                writePixelUncompressed(lzwP, colorIndex);
        }
        ++nRowsDone;
    }
    /* Gif is no good with no pixels; fortunately, that's impossible: */
    assert(nRowsDone > 0);

    lzw_flush(lzwP);

    pnm_freepamrow(tuplerow);

    lzw_destroy(lzwP);
}



static void
writeGlobalColorMap(FILE *              const ofP,
                    const struct Cmap * const cmapP,
                    unsigned int        const bitsPerPixel) {
/*----------------------------------------------------------------------------
  Write out the Global Color Map

  Note that the Global Color Map is always a power of two colors
  in size, but *cmapP could be smaller than that.  So we pad with
  black.
-----------------------------------------------------------------------------*/
    unsigned int const colorMapSize = 1 << bitsPerPixel;

    struct pam pam;
    unsigned int i;
    tuple tupleRgb255;

    if (verbose)
        pm_message("Writing %u-entry global colormap for %u colors",
                   colorMapSize, cmapP->cmapSize);

    pam = cmapP->pam;
    pam.size = PAM_STRUCT_SIZE(allocation_depth);
    pam.len = pam.size;
    pnm_setminallocationdepth(&pam, 3);

    tupleRgb255 = pnm_allocpamtuple(&pam);

    for (i = 0; i < colorMapSize; ++i) {
        if (i < cmapP->cmapSize) {
            tuple const color = cmapP->color[i];

            assert(i < cmapP->cmapSize);

            pnm_scaletuple(&pam, tupleRgb255, color, 255);
            pnm_maketuplergb(&pam, tupleRgb255);

            fputc(tupleRgb255[PAM_RED_PLANE], ofP);
            fputc(tupleRgb255[PAM_GRN_PLANE], ofP);
            fputc(tupleRgb255[PAM_BLU_PLANE], ofP);
        } else {
            fputc(0, ofP);
            fputc(0, ofP);
            fputc(0, ofP);
        }
    }
    pnm_freepamtuple(tupleRgb255);
}



static void
writeGifHeader(FILE *              const ofP,
               unsigned int        const width,
               unsigned int        const height,
               unsigned int        const background,
               unsigned int        const bitsPerPixel,
               const struct Cmap * const cmapP,
               char                const comment[],
               float               const aspect) {

    unsigned int const resolution = bitsPerPixel;

    unsigned char b;

    /* Write the Magic header */
    if (cmapP->haveTransparent || comment || aspect != 1.0 )
        fwrite("GIF89a", 1, 6, ofP);
    else
        fwrite("GIF87a", 1, 6, ofP);

    /* Write out the screen width and height */
    Putword(width,  ofP);
    Putword(height, ofP);

    /* Indicate that there is a global color map */
    b = 0x80;       /* Yes, there is a color map */

    /* OR in the resolution */
    b |= (resolution - 1) << 4;

    /* OR in the Bits per Pixel */
    b |= (bitsPerPixel - 1);

    /* Write it out */
    fputc(b, ofP);

    /* Write out the Background color */
    assert((unsigned char)background == background);
    fputc(background, ofP);

    {
        int const aspectValue = aspect == 1.0 ? 0 : ROUND(aspect * 64) - 15;
        assert(0 <= aspectValue && aspectValue <= 255);
        fputc(aspectValue, ofP);
    }
    writeGlobalColorMap(ofP, cmapP, bitsPerPixel);

    if (cmapP->haveTransparent)
        writeTransparentColorIndexExtension(ofP, cmapP->transparent);

    if (comment)
        writeCommentExtension(ofP, comment);
}



static void
writeImageHeader(FILE *       const ofP,
                 unsigned int const leftOffset,
                 unsigned int const topOffset,
                 unsigned int const gWidth,
                 unsigned int const gHeight,
                 bool         const gInterlace,
                 unsigned int const initCodeSize) {

    Putword(leftOffset, ofP);
    Putword(topOffset,  ofP);
    Putword(gWidth,     ofP);
    Putword(gHeight,    ofP);

    /* Write out whether or not the image is interlaced */
    if (gInterlace)
        fputc(0x40, ofP);
    else
        fputc(0x00, ofP);

    /* Write out the initial code size */
    fputc(initCodeSize, ofP);
}



static void
reportImageInfo(bool         const interlace,
                unsigned int const background,
                unsigned int const bitsPerPixel) {

    if (verbose) {
        if (interlace)
            pm_message("interlaced");
        else
            pm_message("not interlaced");
        pm_message("Background color index = %u", background);
        pm_message("%u bits per pixel", bitsPerPixel);
    }
}



static void
gifEncode(struct pam *  const pamP,
          FILE *        const ofP,
          pm_filepos    const rasterPos,
          bool          const gInterlace,
          int           const background,
          unsigned int  const bitsPerPixel,
          struct Cmap * const cmapP,
          char          const comment[],
          float         const aspect,
          bool          const lzw) {

    unsigned int const leftOffset = 0;
    unsigned int const topOffset  = 0;

    unsigned int const initCodeSize = bitsPerPixel <= 1 ? 2 : bitsPerPixel;
        /* The initial code size */

    sample const alphaThreshold = (pamP->maxval + 1) / 2;
        /* Levels below this in the alpha plane indicate transparent
           pixels in the output image.
        */

    unsigned int const alphaPlane = pamAlphaPlane(pamP);

    RowReader * rowReaderP;

    reportImageInfo(gInterlace, background, bitsPerPixel);

    if (pamP->width > 65535)
        pm_error("Image width %u too large for GIF format.  (Max 65535)",
                 pamP->width);

    if (pamP->height > 65535)
        pm_error("Image height %u too large for GIF format.  (Max 65535)",
                 pamP->height);

    writeGifHeader(ofP, pamP->width, pamP->height, background,
                   bitsPerPixel, cmapP, comment, aspect);

    /* Write an Image separator */
    fputc(',', ofP);

    writeImageHeader(ofP, leftOffset, topOffset, pamP->width, pamP->height,
                     gInterlace, initCodeSize);

    rowReaderP = rowReader_create(pamP, rasterPos, gInterlace);

    /* Write the actual raster */

    writeRaster(pamP, rowReaderP, alphaPlane, alphaThreshold,
                cmapP, initCodeSize + 1, ofP, lzw);

    rowReader_destroy(rowReaderP);

    /* Write out a zero length data block (to end the series) */
    fputc(0, ofP);

    /* Write the GIF file terminator */
    fputc(';', ofP);
}



static void
reportTransparent(struct Cmap * const cmapP) {

    if (verbose) {
        if (cmapP->haveTransparent) {
            tuple const color = cmapP->color[cmapP->transparent];
            pm_message("Color %u (%lu, %lu, %lu) is transparent",
                       cmapP->transparent,
                       color[PAM_RED_PLANE],
                       color[PAM_GRN_PLANE],
                       color[PAM_BLU_PLANE]);
        } else
            pm_message("No transparent color");
    }
}



static void
computeTransparent(char          const colorarg[],
                   bool          const usingFakeTrans,
                   unsigned int  const fakeTransparent,
                   struct Cmap * const cmapP) {
/*----------------------------------------------------------------------------
   Figure out the color index (index into the colormap) of the color
   that is to be transparent in the GIF.

   colorarg[] is the string that specifies the color the user wants to
   be transparent (e.g. "red", "#fefefe").  Its maxval is the maxval
   of the colormap.  'cmap' is the full colormap except that its
   'transparent' component isn't valid.

   colorarg[] is a standard Netpbm color specification, except that
   may have a "=" prefix, which means it specifies a particular exact
   color, as opposed to without the "=", which means "the color that
   is closest to this and actually in the image."

   colorarg[] null means the color didn't ask for a particular color
   to be transparent.

   Establish no transparent color if colorarg[] specifies an exact
   color and that color is not in the image.  Also issue an
   informational message.

   'usingFakeTrans' means pixels will be transparent because of something
   other than their foreground color, and 'fakeTransparent' is the
   color map index for transparent colors.
-----------------------------------------------------------------------------*/
    if (colorarg) {
        const char * colorspec;
        bool exact;
        tuple transcolor;
        int found;
        int colorindex;

        if (colorarg[0] == '=') {
            colorspec = &colorarg[1];
            exact = TRUE;
        } else {
            colorspec = colorarg;
            exact = FALSE;
        }

        transcolor = pnm_parsecolor(colorspec, cmapP->pam.maxval);
        pnm_lookuptuple(&cmapP->pam, cmapP->tuplehash, transcolor, &found,
                        &colorindex);

        if (found) {
            cmapP->haveTransparent = TRUE;
            cmapP->transparent = colorindex;
        } else if (!exact) {
            cmapP->haveTransparent = TRUE;
            cmapP->transparent = closestColor(transcolor, &cmapP->pam, cmapP);
        } else {
            cmapP->haveTransparent = FALSE;
            pm_message("Warning: specified transparent color "
                       "does not occur in image.");
        }
    } else if (usingFakeTrans) {
        cmapP->haveTransparent = TRUE;
        cmapP->transparent = fakeTransparent;
    } else
        cmapP->haveTransparent = FALSE;

    reportTransparent(cmapP);
}



static unsigned int
sortOrderColor(tuple const tuple) {

    return ((tuple[PAM_RED_PLANE] * MAXCMAPSIZE) +
            tuple[PAM_GRN_PLANE]) * MAXCMAPSIZE +
           tuple[PAM_BLU_PLANE];
}



#ifndef LITERAL_FN_DEF_MATCH
static qsort_comparison_fn sortCompareColor;
#endif

static int
sortCompareColor(const void * const entry1P,
                 const void * const entry2P) {

    struct tupleint * const * const tupleint1PP = entry1P;
    struct tupleint * const * const tupleint2PP = entry2P;

    return (sortOrderColor((*tupleint1PP)->tuple)
            - sortOrderColor((*tupleint2PP)->tuple));
}



#ifndef LITERAL_FN_DEF_MATCH
static qsort_comparison_fn sortCompareGray;
#endif

static int
sortCompareGray(const void * const entry1P,
                const void * const entry2P){

    struct tupleint * const * const tupleint1PP = entry1P;
    struct tupleint * const * const tupleint2PP = entry2P;

    return ((*tupleint1PP)->tuple[0] - (*tupleint2PP)->tuple[0]);
}



static void
sortTupletable(struct pam * const mapPamP,
               unsigned int const colors,
               tupletable   const tuplefreq) {
/*----------------------------------------------------------------------------
   Sort the colormap *cmapP.

   Sort the colormap by red intensity, then by green intensity,
   then by blue intensity.
-----------------------------------------------------------------------------*/

    pm_message("sorting colormap");

    if (mapPamP->depth < 3)
        qsort(tuplefreq, colors, sizeof(tuplefreq[0]), sortCompareGray);
    else
        qsort(tuplefreq, colors, sizeof(tuplefreq[0]), sortCompareColor);

}



static void
addToColormap(struct Cmap *  const cmapP,
              const char *   const colorspec,
              unsigned int * const newIndexP) {
/*----------------------------------------------------------------------------
  Add a new entry to the colormap.  Make the color that specified by
  'colorspec', and return the index of the new entry as *newIndexP.

  'colorspec' is a color specification given by the user, e.g.
  "red" or "rgb:ff/03/0d".  The maxval for this color specification is
  that for the colormap *cmapP.
-----------------------------------------------------------------------------*/
    tuple const transcolor = pnm_parsecolor(colorspec, cmapP->pam.maxval);

    unsigned int const colorIndex = cmapP->cmapSize++;

    cmapP->color[colorIndex] = pnm_allocpamtuple(&cmapP->pam);

    if (cmapP->pam.depth < 3) {
        if (!pnm_rgbtupleisgray(transcolor))
            pm_error("Image is grayscale, but color '%s' is not gray.  "
                     "It is (%lu, %lu, %lu)",
                     colorspec,
                     transcolor[PAM_RED_PLANE],
                     transcolor[PAM_GRN_PLANE],
                     transcolor[PAM_BLU_PLANE]);
        else
            cmapP->color[colorIndex][0] = transcolor[0];
    } else {
        pnm_assigntuple(&cmapP->pam, cmapP->color[colorIndex], transcolor);
    }
    *newIndexP = colorIndex;
}



static void
colormapFromFile(char               const filespec[],
                 unsigned int       const maxcolors,
                 tupletable *       const tupletableP,
                 struct pam *       const mapPamP,
                 unsigned int *     const colorCountP) {
/*----------------------------------------------------------------------------
   Read a colormap from the Netpbm file filespec[].  Return a
   tupletable of the colors in it (which is practically a colormap) as
   *tupletableP and the format of those tuples as *mapPamP.  Return
   the number of colors as *colorsCountP.
-----------------------------------------------------------------------------*/
    FILE * mapfileP;
    tuple ** colors;
    unsigned int colorCount;

    mapfileP = pm_openr(filespec);
    colors = pnm_readpam(mapfileP, mapPamP, PAM_STRUCT_SIZE(tuple_type));
    pm_close(mapfileP);

    pm_message("computing other colormap ...");

    *tupletableP =
        pnm_computetuplefreqtable(mapPamP, colors, maxcolors, &colorCount);

    *colorCountP = colorCount;

    pnm_freepamarray(colors, mapPamP);
}



static void
readAndValidateColormapFromFile(char           const filename[],
                                unsigned int   const maxcolors,
                                tupletable *   const tuplefreqP,
                                struct pam *   const mapPamP,
                                unsigned int * const colorCountP,
                                unsigned int   const nInputComp,
                                sample         const inputMaxval) {
/*----------------------------------------------------------------------------
   Read the colormap from a separate colormap file named filename[],
   and make sure it's consistent with an image with 'nInputComp'
   color components (e.g. 3 for RGB) and a maxval of 'inputMaxval'.
-----------------------------------------------------------------------------*/
    colormapFromFile(filename, maxcolors, tuplefreqP, mapPamP, colorCountP);

    if (mapPamP->depth != nInputComp)
        pm_error("Depth of map file (%u) does not match number of "
                 "color components in input file (%u)",
                 mapPamP->depth, nInputComp);
    if (mapPamP->maxval != inputMaxval)
        pm_error("Maxval of map file (%lu) does not match maxval of "
                 "input file (%lu)", mapPamP->maxval, inputMaxval);
}



static void
computeColormapBw(struct pam *   const pamP,
                  struct pam *   const mapPamP,
                  unsigned int * const colorCountP,
                  tupletable   * const tuplefreqP) {
/*----------------------------------------------------------------------------
  Shortcut for black and white (e.g. PBM).  We know that there are
  only two colors.  Users who know that only one color is present in
  the image should specify -sort at the command line.  Example:

   $ pbmmake -w 600 400 | pamtogif -sort > canvas.gif
-----------------------------------------------------------------------------*/
    tupletable const colormap = pnm_alloctupletable(pamP, 2);

    *mapPamP = *pamP;
    mapPamP->depth = 1;

    colormap[0]->value = 1;
    colormap[0]->tuple[0] = PAM_BLACK;
    colormap[1]->value = 1;
    colormap[1]->tuple[0] = PAM_BW_WHITE;

    *tuplefreqP  = colormap;
    *colorCountP = 2;
}



static void
computeColormapFromInput(struct pam *   const pamP,
                         unsigned int   const maxcolors,
                         unsigned int   const nInputComp,
                         struct pam *   const mapPamP,
                         unsigned int * const colorCountP,
                         tupletable *   const tuplefreqP) {

    tupletable tuplefreq;

    pm_message("computing colormap...");

    tuplefreq = pnm_computetuplefreqtable3(
        pamP, NULL, maxcolors, nInputComp, pamP->maxval, colorCountP);

    *mapPamP = *pamP;
    mapPamP->depth = nInputComp;

    *tuplefreqP = tuplefreq;
}



static void
computeLibnetpbmColormap(struct pam *   const pamP,
                         bool           const haveAlpha,
                         const char *   const mapfile,
                         tuple *        const color,
                         tuplehash *    const tuplehashP,
                         struct pam *   const mapPamP,
                         unsigned int * const colorCountP,
                         bool           const sort) {
/*----------------------------------------------------------------------------
   Compute a colormap, libnetpbm style, for the image described by
   'pamP', which is positioned to the raster.

   If 'mapfile' is non-null, Use the colors in that (Netpbm) file for
   the color map instead of the colors in 'pamP'.

   Return the colormap as color[] and *tuplehashP.  Return the format
   of those tuples as *mapPamP.

   The tuples of the color map have a meaningful depth of 1 (grayscale) or 3
   (color) and *mapPamP reflects that.

   While we're at it, count the colors and validate that there aren't
   too many.  Return the count as *colorCountP.  In determining if there are
   too many, allow one slot for a fake transparency color if 'haveAlpha'
   is true.  If there are too many, issue an error message and abort the
   program.

   'sort' means to sort the colormap by red intensity, then by green
   intensity, then by blue intensity, as opposed to arbitrary order.
-----------------------------------------------------------------------------*/
    unsigned int const maxcolors = haveAlpha ? MAXCMAPSIZE - 1 : MAXCMAPSIZE;
        /* The most colors we can tolerate in the image.  If we have
           our own made-up entry in the colormap for transparency, it
           isn't included in this count.
        */
    unsigned int const nInputComp = haveAlpha ? pamP->depth - 1 : pamP->depth;
        /* Number of color components (not alpha) in the input image */

    unsigned int i;
    tupletable tuplefreq;
    unsigned int colorCount;

    if (mapfile)
        readAndValidateColormapFromFile(mapfile, maxcolors, &tuplefreq,
                                        mapPamP, &colorCount,
                                        nInputComp, pamP->maxval);
    else if (nInputComp == 1 && pamP->maxval == 1 && !sort &&
             pamP->height * pamP->width > 1)
        computeColormapBw(pamP, mapPamP, &colorCount, &tuplefreq);
    else
        computeColormapFromInput(pamP, maxcolors, nInputComp,
                                 mapPamP, &colorCount, &tuplefreq);

    if (tuplefreq == NULL)
        pm_error("too many colors - try doing a 'pnmquant %u'", maxcolors);

    pm_message("%u colors found", colorCount);

    if (sort)
        sortTupletable(mapPamP, colorCount, tuplefreq);

    for (i = 0; i < colorCount; ++i) {
        color[i] = pnm_allocpamtuple(mapPamP);
        pnm_assigntuple(mapPamP, color[i], tuplefreq[i]->tuple);
    }

    /* And make a hash table for fast lookup. */
    *tuplehashP =
        pnm_computetupletablehash(mapPamP, tuplefreq, colorCount);

    *colorCountP = colorCount;

    pnm_freetupletable(mapPamP, tuplefreq);
}



static void
destroyCmap(struct Cmap * const cmapP) {

    unsigned int colorIndex;

    for (colorIndex = 0; colorIndex < cmapP->cmapSize; ++colorIndex)
        pnm_freepamtuple(cmapP->color[colorIndex]);

    pnm_destroytuplehash(cmapP->tuplehash);
}



int
main(int argc, char *argv[]) {
    struct CmdlineInfo cmdline;
    FILE * ifP;
    struct pam pam;
    unsigned int bitsPerPixel;
    pm_filepos rasterPos;

    struct Cmap cmap;
        /* The colormap, with all its accessories */
    unsigned int fakeTransparent;
        /* colormap index of the fake transparency color we're using to
           implement the alpha mask.  Undefined if we're not doing an alpha
           mask.
        */

    pnm_init(&argc, argv);

    parseCommandLine(argc, argv, &cmdline);

    verbose = cmdline.verbose;

    ifP = pm_openr_seekable(cmdline.input_filespec);

    pnm_readpaminit(ifP, &pam, PAM_STRUCT_SIZE(tuple_type));

    pm_tell2(ifP, &rasterPos, sizeof(rasterPos));

    computeLibnetpbmColormap(&pam, !!pamAlphaPlane(&pam), cmdline.mapfile,
                             cmap.color, &cmap.tuplehash,
                             &cmap.pam, &cmap.cmapSize, cmdline.sort);

    assert(cmap.pam.maxval == pam.maxval);

    if (pamAlphaPlane(&pam)) {
        /* Add a fake entry to the end of the colormap for transparency.
           Make its color black.
        */
        addToColormap(&cmap, cmdline.alphacolor, &fakeTransparent);
    }

    bitsPerPixel = cmap.cmapSize == 1 ? 1 : nSignificantBits(cmap.cmapSize-1);

    computeTransparent(cmdline.transparent,
                       !!pamAlphaPlane(&pam), fakeTransparent, &cmap);

    /* All set, let's do it. */
    gifEncode(&pam, stdout, rasterPos,
              cmdline.interlace, 0, bitsPerPixel, &cmap, cmdline.comment,
              cmdline.aspect, !cmdline.nolzw);

    destroyCmap(&cmap);

    pm_close(ifP);
    pm_close(stdout);

    return 0;
}



/*============================================================================
  Original version, named 'ppmgif' was by Jef Poskanzer in 1989, based
  on GIFENCOD by David Rowley <mgardi@watdscu.waterloo.edu>.A Lempel-Zim
  compression based on "compress".

  Switched to use libnetpbm PAM facilities (ergo process PAM images)
  and renamed 'pamtogif' by Bryan Henderson November 2006.

  The non-LZW GIF generation stuff was adapted from the Independent
  JPEG Group's djpeg on 2001.09.29.  In 2006.12 the output subroutines
  were rewritten; now no uncompressed output subroutines are derived from
  the Independent JPEG Group's source code.

  2007.01  Changed sort routine to qsort.  (afu)
  2007.03  Implemented variable hash table size, PBM color table
           shortcut and "-aspect" command line option.   (afu)


  Copyright (C) 1989 by Jef Poskanzer.

  Permission to use, copy, modify, and distribute this software and its
  documentation for any purpose and without fee is hereby granted, provided
  that the above copyright notice appear in all copies and that both that
  copyright notice and this permission notice appear in supporting
  documentation.  This software is provided "as is" without express or
  implied warranty.

  The Graphics Interchange Format(c) is the Copyright property of
  CompuServe Incorporated.  GIF(sm) is a Service Mark property of
  CompuServe Incorporated.
============================================================================*/

