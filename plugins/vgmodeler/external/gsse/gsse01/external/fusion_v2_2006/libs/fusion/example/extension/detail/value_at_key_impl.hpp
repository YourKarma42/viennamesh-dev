/*=============================================================================
    Copyright (c) 2001-2006 Joel de Guzman
    Copyright (c) 2005-2006 Dan Marsden

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(BOOST_FUSION_VALUE_AT_KEY_IMPL_20060223_2025)
#define BOOST_FUSION_VALUE_AT_KEY_IMPL_20060223_2025

namespace fields
{
    struct name;
    struct age;
}

namespace boost { namespace fusion {

    struct example_sequence_tag;

    namespace extension
    {
        template<typename Tag>
        struct value_at_key_impl;

        template<>
        struct value_at_key_impl<example_sequence_tag>
        {
            template<typename Sequence, typename N>
            struct apply;

            template<typename Sequence>
            struct apply<Sequence, fields::name>
            {
                typedef std::string type;
            };

            template<typename Sequence>
            struct apply<Sequence, fields::age>
            {
                typedef int type;
            };
        };
    }
}}

#endif
