/*=============================================================================
    Copyright (c) 2001-2006 Joel de Guzman
    Copyright (c) 2006 Dan Marsden

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(BOOST_FUSION_AT_KEY_20060304_1755)
#define BOOST_FUSION_AT_KEY_20060304_1755

#include <boost/type_traits/is_const.hpp>
#include <boost/fusion/support/tag_of.hpp>
#include <boost/fusion/support/detail/access.hpp>

namespace boost { namespace fusion
{
    namespace extension
    {
        template <typename Tag>
        struct at_key_impl
        {
            template <typename Sequence, typename Key>
            struct apply;
        };
    }

    namespace result_of
    {
        template <typename Sequence, typename Key>
        struct at_key
        {
            typedef typename
                extension::at_key_impl<typename traits::tag_of<Sequence>::type>::
                    template apply<Sequence, Key>::type
            type;
        };
    }

    template <typename Key, typename Sequence>
    inline typename 
        lazy_disable_if<
            is_const<Sequence>
          , result_of::at_key<Sequence, Key>
        >::type
    at_key(Sequence& seq)
    {
        typedef result_of::at_key<Sequence, Key> at_meta;
        return extension::at_key_impl<typename traits::tag_of<Sequence>::type>::
            template apply<Sequence, Key>::call(seq);
    }

    template <typename Key, typename Sequence>
    inline typename result_of::at_key<Sequence const, Key>::type
    at_key(Sequence const& seq)
    {
        typedef result_of::at_key<Sequence const, Key> at_meta;
        return extension::at_key_impl<typename traits::tag_of<Sequence>::type>::
            template apply<Sequence const, Key>::call(seq);
    }
}}

#endif
