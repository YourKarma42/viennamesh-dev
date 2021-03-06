/*=============================================================================
    Copyright (c) 2001-2006 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_DISTANCE_09172005_0721)
#define FUSION_DISTANCE_09172005_0721

#include <boost/fusion/iterator/detail/distance.hpp>
#include <boost/fusion/support/category_of.hpp>

#include <boost/mpl/int.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

namespace boost { namespace fusion
{
    struct random_access_traversal_tag;

    namespace extension
    {
        template <typename Tag>
        struct distance_impl
        {
            // default implementation
            template <typename First, typename Last>
            struct apply : distance_detail::linear_distance<First, Last> 
            {
                typedef typename traits::category_of<First>::type first_category;
                typedef typename traits::category_of<Last>::type last_category;
                BOOST_MPL_ASSERT((is_same<first_category, last_category>));
                BOOST_MPL_ASSERT_NOT((is_same<first_category, random_access_traversal_tag>));
            };
        };
    }

    namespace result_of
    {
        template <typename First, typename Last>
        struct distance
            : extension::distance_impl<typename First::ftag>:: template apply<First, Last>
        {};
    }
        
    template <typename First, typename Last>
    inline typename result_of::distance<First, Last>::type
    distance(First const& a, Last const& b)
    {
        return extension::distance_impl<typename First::ftag>::
            template apply<First, Last>::call(a, b);
    }
}}

#endif
