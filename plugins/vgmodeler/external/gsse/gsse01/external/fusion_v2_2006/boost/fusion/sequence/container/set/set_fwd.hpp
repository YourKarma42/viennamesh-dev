/*=============================================================================
    Copyright (c) 2005 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_SET_FORWARD_09162005_1102)
#define FUSION_SET_FORWARD_09162005_1102

#include <boost/fusion/sequence/container/set/limits.hpp>
#include <boost/preprocessor/repetition/enum_params_with_a_default.hpp>

namespace boost { namespace fusion
{
    struct void_;

    template <
        BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT(
            FUSION_MAX_SET_SIZE, typename T, void_)
    >
    struct set;
}}

#endif
