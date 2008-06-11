/* pbmtext.c - render text into a bitmap
**
** Copyright (C) 1991 by Jef Poskanzer.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty.
*/

#define _BSD_SOURCE 1      /* Make sure strdup() is in string.h */
#define _XOPEN_SOURCE 500  /* Make sure strdup() is in string.h */

#include <string.h>
#include <limits.h>
#include <assert.h>

#include "mallocvar.h"
#include "shhopt.h"
#include "pbm.h"
#include "pbmfont.h"

struct cmdlineInfo {
    /* All the information the user supplied in the command line,
       in a form easy for the program to use.
    */
    const char * text;    /* text from command line or NULL if none */
    const char * font;    /* -font option value or NULL if none */
    const char * builtin; /* -builtin option value or NULL if none */
    unsigned int dump;   
        /* undocumented dump option for installing a new built-in font */
    float space;   /* -space option value or default */
    unsigned int width;     /* -width option value or zero */
    int lspace;    /* lspace option value or default */
    unsigned int nomargins;     /* -nomargins */
    unsigned int verbose;
};




static void
parseCommandLine(int argc, char ** argv,
                 struct cmdlineInfo * const cmdlineP) {
/*----------------------------------------------------------------------------
   Note that the file spec array we return is stored in the storage that
   was passed to us as the argv array.
-----------------------------------------------------------------------------*/
    optEntry *option_def = malloc(100*sizeof(optEntry));
        /* Instructions to OptParseOptions3 on how to parse our options.
         */
    optStruct3 opt;

    unsigned int option_def_index;

    option_def_index = 0;   /* incremented by OPTENTRY */
    OPTENT3(0, "font",      OPT_STRING, &cmdlineP->font, NULL,        0);
    OPTENT3(0, "builtin",   OPT_STRING, &cmdlineP->builtin, NULL,     0);
    OPTENT3(0, "dump",      OPT_FLAG,   NULL, &cmdlineP->dump,        0);
    OPTENT3(0, "space",     OPT_FLOAT,  &cmdlineP->space, NULL,       0);
    OPTENT3(0, "width",     OPT_UINT,   &cmdlineP->width, NULL,       0);
    OPTENT3(0, "lspace",    OPT_INT,    &cmdlineP->lspace, NULL,      0);
    OPTENT3(0, "nomargins", OPT_FLAG,   NULL, &cmdlineP->nomargins,   0);
    OPTENT3(0, "verbose",   OPT_FLAG,   NULL, &cmdlineP->verbose,     0);

    /* Set the defaults */
    cmdlineP->font = NULL;
    cmdlineP->builtin = NULL;
    cmdlineP->space = 0.0;
    cmdlineP->width = 0;
    cmdlineP->lspace = 0;

    opt.opt_table = option_def;
    opt.short_allowed = FALSE;  /* We have no short (old-fashioned) options */
    opt.allowNegNum = FALSE;  /* We have no parms that are negative numbers */

    optParseOptions3(&argc, argv, opt, sizeof(opt), 0);
        /* Uses and sets argc, argv, and some of *cmdlineP and others. */

    if (argc-1 == 0)
        cmdlineP->text = NULL;
    else {
        char *text;
        int i;
        int totaltextsize;

        totaltextsize = 1;  /* initial value */
        
        text = malloc(totaltextsize);  /* initial allocation */
        text[0] = '\0';
        
        for (i = 1; i < argc; ++i) {
            if (i > 1) {
                totaltextsize += 1;
                text = realloc(text, totaltextsize);
                if (text == NULL)
                    pm_error("out of memory allocating space for input text");
                strcat(text, " ");
            } 
            totaltextsize += strlen(argv[i]);
            text = realloc(text, totaltextsize);
            if (text == NULL)
                pm_error("out of memory allocating space for input text");
            strcat(text, argv[i]);
        }
        cmdlineP->text = text;
    }
}



static void
reportFont(struct font * const fontP) {

    unsigned int n;
    unsigned int c;

    pm_message("FONT:");
    pm_message("  character dimensions: %uw x %uh",
               fontP->maxwidth, fontP->maxheight);
    pm_message("  Additional vert white space: %d pixels", fontP->y);

    for (c = 0, n = 0; c < ARRAY_SIZE(fontP->glyph); ++c)
        if (fontP->glyph[c])
            ++n;

    pm_message("  # characters: %u", n);
}



static void
computeFont(struct cmdlineInfo const cmdline,
            struct font **     const fontPP) {

    struct font * fontP;

    if (cmdline.font)
        fontP = pbm_loadfont(cmdline.font);
    else {
        if (cmdline.builtin)
            fontP = pbm_defaultfont(cmdline.builtin);
        else
            fontP = pbm_defaultfont("bdf");
    }

    if (cmdline.verbose)
        reportFont(fontP);

    if (cmdline.dump) {
        pbm_dumpfont(fontP);
        exit(0);
    }
    *fontPP = fontP;
}



struct text {
    char **      textArray;  /* malloc'ed */
    unsigned int allocatedLineCount;
    unsigned int lineCount;
};



static void
allocTextArray(struct text * const textP,
               unsigned int  const maxLineCount,
               unsigned int  const maxColumnCount) {

    unsigned int line;

    textP->allocatedLineCount = maxLineCount;

    MALLOCARRAY_NOFAIL(textP->textArray, maxLineCount);
    
    for (line = 0; line < maxLineCount; ++line)
        MALLOCARRAY_NOFAIL(textP->textArray[line], maxColumnCount+1);

    textP->lineCount = 0;
}



static void
freeTextArray(struct text const text) {

    unsigned int line;

    for (line = 0; line < text.allocatedLineCount; ++line)
        free((char **)text.textArray[line]);

    free(text.textArray);
}



static void
fixControlChars(const char *  const input,
                struct font * const fontP,
                const char ** const outputP) {
/*----------------------------------------------------------------------------
   Return a translation of input[] that can be rendered as glyphs in
   the font 'fontP'.  Return it as newly malloced *outputP.

   Expand tabs to spaces.

   Remove any trailing newline.  (But leave intermediate ones as line
   delimiters).

   Turn anything that isn't a code point in the font to a single space
   (which isn't guaranteed to be in the font either, of course).
-----------------------------------------------------------------------------*/
    /* We don't know in advance how big the output will be because of the
       tab expansions.  So we make sure before processing each input
       character that there is space in the output buffer for a worst
       case tab expansion, plus a terminating NUL, reallocating as
       necessary.  And we originally allocate enough for the entire line
       assuming no tabs.
    */

    unsigned int const tabSize = 8;

    unsigned int inCursor, outCursor;
    char * output;      /* Output buffer.  Malloced */
    size_t outputSize;  /* Currently allocated size of 'output' */

    outputSize = strlen(input) + 1 + tabSize;
        /* Leave room for one worst case tab expansion and NUL terminator */
    MALLOCARRAY(output, outputSize);

    if (output == NULL)
        pm_error("Couldn't allocate %u bytes for a line of text.", outputSize);

    for (inCursor = 0, outCursor = 0; input[inCursor] != '\0'; ++inCursor) {
        if (outCursor + 1 + tabSize > outputSize) {
            outputSize = outCursor + 1 + 4 * tabSize;
            REALLOCARRAY(output, outputSize);
            if (output == NULL)
                pm_error("Couldn't allocate %u bytes for a line of text.",
                         outputSize);
        }
        if (input[inCursor] == '\n' && input[inCursor+1] == '\0') {
            /* This is a terminating newline.  We don't do those. */
        } else if (input[inCursor] == '\t') { 
            /* Expand this tab into the right number of spaces. */
            unsigned int const nextTabStop =
                (outCursor + tabSize) / tabSize * tabSize;

            while (outCursor < nextTabStop)
                output[outCursor++] = ' ';
        } else if (!fontP->glyph[(unsigned char)input[inCursor]]) {
            /* Turn this unknown char into a single space. */
            output[outCursor++] = ' ';
        } else
            output[outCursor++] = input[inCursor];

        assert(outCursor <= outputSize);
    }
    output[outCursor++] = '\0';

    assert(outCursor <= outputSize);

    *outputP = output;
}



static void
fill_rect(bit** const bits, 
          int   const row0, 
          int   const col0, 
          int   const height, 
          int   const width, 
          bit   const color) {

    int row;

    for (row = row0; row < row0 + height; ++row) {
        int col;
        for (col = col0; col < col0 + width; ++col)
            bits[row][col] = color;
    }
}



static void
get_line_dimensions(const char line[], const struct font * const font_p, 
                    const float intercharacter_space,
                    double * const bwidP, int * const backup_space_needed_p) {
/*----------------------------------------------------------------------------
   Determine the width in pixels of the line of text line[] in the font
   *font_p, and return it as *bwidP.  Also determine how much of this
   width goes to the left of the nominal starting point of the line because
   the first character in the line has a "backup" distance.  Return that
   as *backup_space_needed_p.
-----------------------------------------------------------------------------*/
    int cursor;  /* cursor into the line of text */
    double accumulatedIcs;
        /* accumulated intercharacter space so far in the line we are 
           stepping through.  Because the intercharacter space might not be
           an integer, we accumulate it here and realize full pixels whenever
           we have more than one pixel.  Note that this can be negative
           (which means were crowding, rather than spreading, text).
        */
    double bwid;
    bool no_chars_yet; 
        /* We haven't seen any renderable characters yet in the line. */
    struct glyph * lastGlyphP;
        /* Glyph of last character processed so far.  Undefined if
           'no_chars_yet'.
        */

    no_chars_yet = TRUE;   /* initial value */
    accumulatedIcs = 0.0;  /* initial value */
    bwid = 0.0;  /* initial value */
    
    for (cursor = 0; line[cursor] != '\0'; cursor++) {
        struct glyph * const glyphP = 
            font_p->glyph[(unsigned char)line[cursor]];

        if (glyphP) {
            if (no_chars_yet) {
                no_chars_yet = FALSE;
                if (glyphP->x < 0) 
                    *backup_space_needed_p = -glyphP->x;
                else {
                    *backup_space_needed_p = 0;
                    bwid += glyphP->x;
                }
            } else {
                /* handle extra intercharacter space (-space option) */
                accumulatedIcs += intercharacter_space;
                if (accumulatedIcs >= INT_MAX)
                    pm_error("Image width too large.");
                if (accumulatedIcs <= INT_MIN)
                    pm_error("Absurdly large negative -space value.");
                {
                    int const fullPixels = (int) accumulatedIcs;
                    bwid           += fullPixels;
                    accumulatedIcs -= fullPixels;
                }
            }
            lastGlyphP = glyphP;
            bwid += glyphP->xadd;
        }
    }
    if (no_chars_yet)
        /* Line has no renderable characters */
        *backup_space_needed_p = 0;
    else {
        /* Line has at least one renderable character.
           Recalculate width of last character in line so it ends
           right at the right edge of the glyph (no extra space to
           anticipate another character).
        */
        bwid -= lastGlyphP->xadd;
        bwid += lastGlyphP->width + lastGlyphP->x;
    }
    if (bwid > INT_MAX)
        pm_error("Image width too large.");
    else
        *bwidP = bwid; 
}



static void
insert_character(const struct glyph * const glyph, 
                 int                  const toprow, 
                 int                  const leftcol,
                 bit **               const bits) {
/*----------------------------------------------------------------------------
   Insert one character (whose glyph is 'glyph') into the image bits[].
   Its top left corner shall be row 'toprow', column 'leftcol'.
-----------------------------------------------------------------------------*/

    int glyph_y, glyph_x;  /* position within the glyph */

    for (glyph_y = 0; glyph_y < glyph->height; glyph_y++) {
        for (glyph_x = 0; glyph_x < glyph->width; glyph_x++) {
            if (glyph->bmap[glyph_y * glyph->width + glyph_x])
                bits[toprow+glyph_y][leftcol+glyph->x+glyph_x] = 
                    PBM_BLACK;
        }
    }
}    



static void
insert_characters(bit **        const bits, 
                  struct text   const lp,
                  struct font * const fontP, 
                  int           const topmargin, 
                  int           const leftmargin,
                  float         const intercharacter_space,
                  int           const lspace) {
/*----------------------------------------------------------------------------
   Render the text 'lp' into the image 'bits' using font *fontP and
   putting 'intercharacter_space' pixels between characters and
   'lspace' pixels between the lines.
-----------------------------------------------------------------------------*/
    int line;  /* Line number in input text */

    for (line = 0; line < lp.lineCount; ++line) {
        int row;  /* row in image of top of current typeline */
        int leftcol;  /* Column in image of left edge of current glyph */
        int cursor;  /* cursor into a line of input text */
        float accumulated_ics;
            /* accumulated intercharacter space so far in the line we
               are building.  Because the intercharacter space might
               not be an integer, we accumulate it here and realize
               full pixels whenever we have more than one pixel. 
            */

        row = topmargin + line * (fontP->maxheight + lspace);
        leftcol = leftmargin;
        accumulated_ics = 0.0;  /* initial value */
    
        for (cursor = 0; lp.textArray[line][cursor] != '\0'; ++cursor) {
            unsigned int const glyphIndex = 
                (unsigned char)lp.textArray[line][cursor];
            struct glyph* glyph;   /* the glyph for this character */

            glyph = fontP->glyph[glyphIndex];
            if (glyph != NULL) {
                const int toprow = row + fontP->maxheight + fontP->y 
                    - glyph->height - glyph->y;
                    /* row number in image of top row in glyph */
                
                insert_character(glyph, toprow, leftcol, bits);

                leftcol += glyph->xadd;
                {
                    /* handle extra intercharacter space (-space option) */
                    int full_pixels;  /* integer part of accumulated_ics */
                    accumulated_ics += intercharacter_space;
                    full_pixels = (int) accumulated_ics;
                    if (full_pixels > 0) {
                        leftcol += full_pixels;
                        accumulated_ics -= full_pixels;
                    }
                }
            }
        }
    }
}


struct outputTextCursor {
    struct text text;
        /* The output text.  The lineCount field of this represents
           the number of lines we have completed.  The line after that
           is the one we are currently filling.
        */
    unsigned int maxWidth;
        /* A line of output can't be wider than this many pixels */
    float intercharacterSpace;
        /* The amount of extra space, in characters, that should be added
           between every two characters (Pbmtext -space option)
        */
    unsigned int columnNo;
        /* The column Number (starting at 0) in the current line that we are
           filling where the next character goes.
        */
    bool noCharsYet;
        /* We haven't put any renderable characters yet in the
           output line. 
        */
    unsigned int widthSoFar;
        /* The accumulated width, in pixels, of all the characters now
           in the current output line 
        */
    float accumulatedIcs;
        /* accumulated intercharacter space so far in the line we
           are stepping through.  Because the intercharacter space
           might not be an integer, we accumulate it here and
           realize full pixels whenever we have more than one
           pixel.  Note that this is negative if we're crowding, rather
           than spreading, characters.
        */
};



static void
initializeFlowedOutputLine(struct outputTextCursor * const cursorP) {

    cursorP->columnNo = 0;
    cursorP->noCharsYet = TRUE;
    cursorP->widthSoFar = 0.0;
    cursorP->accumulatedIcs = 0.0;
}



static void
initializeFlowedOutput(struct outputTextCursor * const cursorP,
                       unsigned int              const maxLines,
                       unsigned int              const maxWidth,
                       float                     const intercharacterSpace) {
    
    allocTextArray(&cursorP->text, maxLines, maxWidth);
    cursorP->maxWidth = maxWidth;
    cursorP->intercharacterSpace = intercharacterSpace;
    initializeFlowedOutputLine(cursorP);
}



static void
finishOutputLine(struct outputTextCursor * const cursorP) {

    if (cursorP->text.lineCount < cursorP->text.allocatedLineCount) {
        char * const currentLine = 
            cursorP->text.textArray[cursorP->text.lineCount];
        currentLine[cursorP->columnNo++] = '\0';
        ++cursorP->text.lineCount;
    }
}



static void
placeCharacterInOutput(char                      const lastch,
                       struct font *             const fontP, 
                       struct outputTextCursor * const cursorP) {
/*----------------------------------------------------------------------------
   Place a character of text in the text array at the position indicated
   by *cursorP, keeping track of what space this character will occupy
   when this text array is ultimately rendered using font *fontP.

   Note that while we compute how much space the character will take when
   rendered, we don't render it.
-----------------------------------------------------------------------------*/
    if (cursorP->text.lineCount < cursorP->text.allocatedLineCount) {
        unsigned int const glyphIndex = (unsigned char)lastch;
        if (fontP->glyph[glyphIndex]) {
            if (cursorP->noCharsYet) {
                cursorP->noCharsYet = FALSE;
                if (fontP->glyph[glyphIndex]->x > 0) 
                    cursorP->widthSoFar += fontP->glyph[glyphIndex]->x;
            } else {
                /* handle extra intercharacter space (-space option) */
                cursorP->accumulatedIcs += cursorP->intercharacterSpace;
                {
                    int const fullPixels = (int)cursorP->accumulatedIcs;
                    cursorP->widthSoFar     += fullPixels;
                    cursorP->accumulatedIcs -= fullPixels;
                }
            }
            cursorP->widthSoFar += fontP->glyph[glyphIndex]->xadd;
        }
        if (cursorP->widthSoFar < cursorP->maxWidth) {
            char * const currentLine = 
                cursorP->text.textArray[cursorP->text.lineCount];
            currentLine[cursorP->columnNo++] = lastch;
        } else {
            /* Line is full; finish it off, start the next one, and
               place the character there.
            */
            /* TODO: We really should back up to the previous white space
               character and move the rest of the line to the next line
            */
            finishOutputLine(cursorP);
            initializeFlowedOutputLine(cursorP);
            placeCharacterInOutput(lastch, fontP, cursorP);
        }
    }
}



static void
flowText(struct text    const inputText,
         int            const width, 
         struct font *  const fontP, 
         float          const intercharacterSpace,
         struct text *  const outputTextP) {
    
    unsigned int const maxLineCount = 50;

    unsigned int inputLine;  
        /* Input line number on which we are currently working */
    struct outputTextCursor outputCursor;

    for (inputLine = 0; inputLine < inputText.lineCount; ++inputLine) {
        unsigned int incursor;   /* cursor into the line we are reading */

        initializeFlowedOutput(&outputCursor, maxLineCount,
                               width, intercharacterSpace);
        
        for (incursor = 0; 
             inputText.textArray[inputLine][incursor] != '\0'; 
             ++incursor)
            placeCharacterInOutput(inputText.textArray[inputLine][incursor],
                                   fontP, &outputCursor);
        finishOutputLine(&outputCursor);
    }
    *outputTextP = outputCursor.text;
}



static void
truncateText(struct text   const inputText, 
             unsigned int  const width, 
             struct font * const fontP, 
             float         const intercharacterSpace,
             struct text * const outputTextP) {

    struct text truncatedText;
    int line;  /* Line number on which we are currently working */

    allocTextArray(&truncatedText, inputText.lineCount, width);

    for (line = 0; line < inputText.lineCount; ++line){
        int cursor;  /* cursor into the line of text */
        unsigned char lastch;  /* line[cursor] */
        int widthSoFar;
            /* How long the line we've built, in pixels, is so far */
        float accumulatedIcs;
        /* accumulated intercharacter space so far in the line we are 
           stepping through.  Because the intercharacter space might not be
           an integer, we accumulate it here and realize full pixels whenever
           we have more than one pixel.  Note that this is negative if we're
           crowding, not spreading, characters.
        */

        int noCharsYet; 
        /* logical: we haven't seen any renderable characters yet in 
           the line.
        */
        noCharsYet = TRUE;   /* initial value */
        widthSoFar = 0;  /* initial value */
        accumulatedIcs = 0.0;  /* initial value */

        truncatedText.textArray[line][0] = '\0';  /* Start with empty line */

        for (cursor = 0; 
             inputText.textArray[line][cursor] != '\0' && widthSoFar < width; 
             cursor++) {
            lastch = inputText.textArray[line][cursor];
            if (fontP->glyph[(unsigned char)lastch]) {
                if (noCharsYet) {
                    noCharsYet = FALSE;
                    if (fontP->glyph[lastch]->x > 0) 
                        widthSoFar += fontP->glyph[lastch]->x;
                } else {
                    /* handle extra intercharacter space (-space option) */
                    accumulatedIcs += intercharacterSpace;
                    {
                        int const fullPixels = (int) intercharacterSpace;
                        widthSoFar     += fullPixels;
                        accumulatedIcs -= fullPixels;
                    }
                }
                widthSoFar += fontP->glyph[lastch]->xadd;
            }
            if (widthSoFar < width) {
                truncatedText.textArray[line][cursor] = 
                    inputText.textArray[line][cursor];
                truncatedText.textArray[line][cursor+1] = '\0';
            }
        }
    }
    truncatedText.lineCount = inputText.lineCount;
    *outputTextP = truncatedText;
}



static void
getText(const char          cmdline_text[], 
        struct font * const fontP,
        struct text * const input_textP) {

    struct text input_text;

    if (cmdline_text) {
        MALLOCARRAY_NOFAIL(input_text.textArray, 1);
        input_text.allocatedLineCount = 1;
        input_text.lineCount = 1;
        fixControlChars(cmdline_text, fontP,
                        (const char**)&input_text.textArray[0]);
    } else {
        /* Read text from stdin. */

        unsigned int maxlines;  
            /* Maximum number of lines for which we presently have space
               in the text array 
            */
        char buf[5000];
        char ** text_array;
        unsigned int lineCount;

        maxlines = 50;  /* initial value */
        MALLOCARRAY_NOFAIL(text_array, maxlines);
        
        lineCount = 0;  /* initial value */
        while (fgets(buf, sizeof(buf), stdin) != NULL) {
            if (strlen(buf) + 1 >= sizeof(buf))
                pm_error("A line of input text is longer than %u characters."
                         "Cannot process.", sizeof(buf)-1);
            if (lineCount >= maxlines) {
                maxlines *= 2;
                REALLOCARRAY(text_array, maxlines);
                if (text_array == NULL)
                    pm_error("out of memory");
            }
            fixControlChars(buf, fontP, (const char **)&text_array[lineCount]);
            if (text_array[lineCount] == NULL)
                pm_error("out of memory");
            ++lineCount;
        }
        input_text.textArray = text_array;
        input_text.lineCount = lineCount;
        input_text.allocatedLineCount = lineCount;
    }
    *input_textP = input_text;
}



static void
computeImageHeight(struct text         const formattedText, 
                   const struct font * const fontP,
                   int                 const interlineSpace,
                   unsigned int        const vmargin,
                   unsigned int *      const rowsP) {

    if (interlineSpace < 0 && fontP->maxheight < -interlineSpace)
        pm_error("-lspace value (%d) negative and exceeds font height.",
                 interlineSpace);     
    else {
        double const rowsD = 2 * (double) vmargin + 
            (double) formattedText.lineCount * fontP->maxheight + 
            (double) (formattedText.lineCount-1) * interlineSpace;
        
        if (rowsD > INT_MAX-10)
            pm_error("Image height too large.");
        else
            *rowsP = (unsigned int) rowsD;
    }
}



static void
computeImageWidth(struct text         const formattedText, 
                  const struct font * const fontP,
                  float               const intercharacterSpace,
                  unsigned int        const hmargin,
                  unsigned int *      const colsP,
                  int *               const maxleftbP) {

    if (intercharacterSpace < 0 && fontP->maxwidth < -intercharacterSpace)
        pm_error("-space value (%f) negative; exceeds font width.",
                 intercharacterSpace);     
    else {
        /* Find the widest line, and the one that backs up the most past
           the nominal start of the line.
        */
    
        unsigned int line;
        double maxwidth;
        int maxleftb;
        double colsD;

        for (line = 0, maxwidth = 0.0, maxleftb = 0;
             line < formattedText.lineCount;
             ++line) {

            double bwid;
            int backupSpaceNeeded;
            
            get_line_dimensions(formattedText.textArray[line], fontP,
                                intercharacterSpace,
                                &bwid, &backupSpaceNeeded);
            
            maxwidth = MAX(maxwidth, bwid);
            maxleftb = MAX(maxleftb, backupSpaceNeeded);
        }
        colsD = 2 * (double) hmargin + (double) maxwidth;
    
        if (colsD > INT_MAX-10)
            pm_error("Image width too large.");
        else
            *colsP = (unsigned int) colsD;
    
        *maxleftbP = maxleftb;
    }
}



int
main(int argc, char *argv[]) {

    struct cmdlineInfo cmdline;
    bit ** bits;
    unsigned int rows, cols;
    struct font * fontP;
    unsigned int vmargin, hmargin;
    struct text inputText;
    struct text formattedText;
    int maxleftb;

    pbm_init(&argc, argv);

    parseCommandLine(argc, argv, &cmdline);
    
    computeFont(cmdline, &fontP);

    getText(cmdline.text, fontP, &inputText);
       
    if (cmdline.nomargins) {
        vmargin = 0;
        hmargin = 0;
    } else {
        if (inputText.lineCount == 1) {
            vmargin = fontP->maxheight / 2;
            hmargin = fontP->maxwidth;
        } else {
            vmargin = fontP->maxheight;
            hmargin = 2 * fontP->maxwidth;
        }
    }
    
    if (cmdline.width > 0) {
        if (cmdline.width > INT_MAX -10)
            pm_error("-width value too large: %u", cmdline.width);
            
        /* Flow or truncate lines to meet user's width request */
        if (inputText.lineCount == 1) 
            flowText(inputText, cmdline.width, fontP, cmdline.space,
                     &formattedText);
        else
            truncateText(inputText, cmdline.width, fontP, cmdline.space,
                         &formattedText);
        freeTextArray(inputText);
    } else
        formattedText = inputText;
        
    if (formattedText.lineCount == 0)
        pm_error("No input text.");
    
    computeImageHeight(formattedText, fontP, cmdline.lspace, vmargin,
                       &rows);

    computeImageWidth(formattedText, fontP, cmdline.space, hmargin,
                      &cols, &maxleftb);

    if (cols == 0 || rows == 0)
        pm_error("Input is all whitespace and/or non-renderable characters.");

    bits = pbm_allocarray(cols, rows);

    /* Fill background with white */
    fill_rect(bits, 0, 0, rows, cols, PBM_WHITE);

    /* Put the text in  */
    insert_characters(bits, formattedText, fontP, vmargin, hmargin + maxleftb, 
                      cmdline.space, cmdline.lspace);

    pbm_writepbm(stdout, bits, cols, rows, 0);

    pbm_freearray(bits, rows);

    freeTextArray(formattedText);
    pm_close(stdout);

    return 0;
}
