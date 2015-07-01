
#ifndef NTL_config__H
#define NTL_config__H

/*************************************************************************

                          NTL Configuration File
                          ----------------------

This file may be modified prior to building NTL so as to customize
how code is generated to improve performance.  You can safely
ignore this file, and NTL will still work correctly.
Only performance is affected.

All NTL include (".h") files include this file.
By setting these flags here, instead of on the compiler command line,
it is easier to guarantee that NTL library and client code use
consistent settings.  


                                How to do it
                                ------------

You override NTL's default code generation strategies by setting
various flags, as described below.  To set a flag, just replace the
pre-processor directive 'if 0' by 'if 1' for that flag, 
which causes the appropriate macro to be defined.  Of course, 
to unset a flag, just replace the 'if 1' by an 'if 0'.

 *************************************************************************/





/*************************************************************************
 *
 * Group 1: macros which override the default implementation of 
 * long integer arithmetic.
 * On some machines, these flags can yield a significant performance gain.
 * Note that you may only specify at most one of the flags in this group.
 *
 *************************************************************************/

#if 0
#define NTL_LONG_LONG

/*
 *   RECOMMENDED FOR Linux/Pentium PLATFORMS
 *
 *   For platforms that support it, this flag can be set to cause
 *   the low-level multiplication code to use the type "long long",
 *   which on some platforms yields a significant performance gain,
 *   but on others, it can yield no improvement and can even
 *   slow things down.
 *   The only platform where I know this helps is Linux/Pentium.
 *   See below (NTL_LONG_LONG_TYPE) for how to use a type name 
 *   other than "long long".
 *   To re-build after changing this flag:  touch lip.c; make ntl.a
 */

#elif 0
#define NTL_AVOID_FLOAT

/*
 *   RECOMMENDED FOR AIX/PowerPC PLATFORMS
 *
 *   On machines with slow floating point or---more comminly---slow int/float
 *   conversions, this flag can lead to faster code.
 *   I get better code on a PowerPC using this flag.
 *   You would think that NTL_LONG_LONG would be better,
 *   but many compilers generate really stupid "long long" code.
 *   If you set NTL_AVOID_FLOAT, you should probably also
 *   set NTL_TBL_REM (see below).
 *   To re-build after changing this flag:  touch lip.c; make ntl.a
 */

#elif 0
#define NTL_SINGLE_MUL 

/*   This was developed originally to improve performance on
 *   ancient Sparc stations that did not have built-in integer mul
 *   instructions.  Unless you have such an old-timer, I would not
 *   recommend using this option.  This option only works on
 *   32-bit machines with IEEE floating point, and is not truly
 *   portable.  If you use this option, you get a 26-bit radix.
 *   To re-build after changing this flag: make clobber; make ntl.a
 */

#endif



/*************************************************************************
 *
 * Group 2: some useful flags affecting code generation
 *
 *************************************************************************/



#if 0
#define NTL_TBL_REM

/*
 *   RECOMMENDED FOR AIX/PowerPC PLATFORMS
 *
 *   With this flag, some divisions are avoided in the
 *   ZZ_pX multiplication routines.  If you use the NTL_AVOID_FLOAT flag,
 *   then you should probably use this one too.
 *   To re-build after changing this flag: touch ZZ_p.c ZZ_pX.c; make ntl.a
 */

#endif

 
#if 0
#define NTL_RANGE_CHECK

/*
 *   This will generate vector subscript range-check code.
 *   Useful for debugging, but it slows things down of course.
 *   To re-build after changing this flag: make clobber; make ntl.a
 */

#endif


/*************************************************************************
 *
 * Group 3: some slightly esoteric things
 *
 *************************************************************************/



#if 0
#define NTL_LONG_LONG_TYPE long long

/*
 *   If you set NTL_LONG_LONG, you may need to override the default
 *   name of this "nonstandard" type.  For example, under MS C++,
 *   the right name is __int64.
 */

#endif


#if 0
#define NTL_CPLUSPLUS_ONLY

/*
 *   It is possible to compile everything using C++ only.
 *   If you want to do this, make CC and CPP in the makefile the same.
 *   You may also want to set this flag, which eliminates some
 *   "C" linkage that is no longer necessary.
 *   However, it should still work without it.
 *   To re-build after changing this flag: make clobber; make ntl.a
 */

#endif


#if 0
#define DNTL_AVOID_BRANCHING

/*
 *   With this option, branches are replaced at several 
 *   key points with equivalent code using shifts and masks.
 *   Recommended for use with RISC architectures, especially
 *   ones with deep pipelines and high branch penalities.
 *   This flag is becoming less helpful as newer machines
 *   have much smaller branch penalties, but still may be worth a try.
 *   To re-build after changing this flag: make clobber; make ntl.a
 */

#endif


#if 0
#define NTL_FFT_PIPELINE

/*
 *   If using NTL_AVOID_BRANCHING, you might want to try this as well.
 *   This causes the FFT routine to use a software pipeline.
 *   To re-build after changing this flag: touch FFT.c; make ntl.a
 */

#endif
 
#if 0
#define NTL_FAST_INT_MUL

/*
 *   Really esoteric.
 *   If using NTL_SINGLE_MUL, and your machine
 *   has a fast integer multiply instruction, this might yield
 *   faster code.  Experiment!
 *   To re-build after changing this flag: make clobber; make ntl.a
 */

#endif


#if 0
#define NTL_NO_INIT_TRANS

/*
 *   Without this flag, NTL uses a special code sequence to avoid
 *   copying large objects in return statements.  However, if your
 *   compiler optimizes away the return of a *named* local object,
 *   this is not necessary, and setting this flag will result
 *   in *slightly* more compact and efficient code.  Although
 *   the emeriging C++ standard allows compilers to perform
 *   this optimization, I know of none that currently do.
 *   Most will avoid copying *temporary* objects in return statements,
 *   and NTL's default code sequence exploits this fact.
 *   To re-build after changing this flag: make clobber; make ntl.a
 */

#endif


#if 0
#define NTL_X86_FIX

/*
 *  Forces the "x86 floating point fix", overriding the default behavior.
 *  By default, NTL will apply the "fix" if it looks like it is
 *  necessary, and if knows how to fix it.
 *  The problem addressed here is that x86 processors sometimes
 *  run in a mode where FP registers have more precision than doubles.
 *  This will cause code in quad_float.c some trouble.
 *  NTL can normally correctly detect the problem, and fix it,
 *  so you shouldn't need to worry about this or the next flag.
 */

#elif 0
#define NTL_NO_X86_FIX
/*
 *  Forces no "x86 floating point fix", overriding the default behavior.
 */

#endif



/*
 *
 *  Version Number
 *
 */

#define NTL_VERSION "3.7a"

#endif
