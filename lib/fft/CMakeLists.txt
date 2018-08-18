cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

set(sources ./r1f2kb.f ./r1f2kf.f ./r1f3kb.f ./r1f3kf.f ./r1f4kb.f ./r1f4kf.f ./r1f5kb.f ./r1f5kf.f
            ./r1fgkb.f ./r1fgkf.f ./rfft1b.f ./rfft1f.f ./rfft1i.f ./rfft2b.f ./rfft2i.f ./rfftb1.f
            ./rfftf1.f ./rffti1.f ./rfftmb.f ./rfftmf.f ./rfftmi.f ./xerfft.f)

add_library(fftpack ${sources})
