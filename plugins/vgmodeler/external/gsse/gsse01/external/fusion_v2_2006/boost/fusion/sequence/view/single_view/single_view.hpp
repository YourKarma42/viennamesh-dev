/*=============================================================================
    Copyright (c) 2001-2006 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_SINGLE_VIEW_05052005_0335)
#define FUSION_SINGLE_VIEW_05052005_0335

#include <boost/fusion/support/detail/access.hpp>
#include <boost/fusion/support/detail/as_fusion_element.hpp>
#include <boost/fusion/support/sequence_base.hpp>
#include <boost/fusion/sequence/view/single_view/single_view_iterator.hpp>
#include <boost/fusion/sequence/view/single_view/detail/begin_impl.hpp>
#include <boost/fusion/sequence/view/single_view/detail/end_impl.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/int.hpp>

namespace boost { namespace fusion
{
    struct single_view_tag;
    struct forward_sequence_tag;
    struct fusion_sequence_tag;

    template <typename T>
    struct single_view : sequence_base<single_view<T> >
    {
        typedef single_view_tag ftag;
        typedef fusion_sequence_tag tag; // this gets picked up by MPL
        typedef forward_sequence_tag category;
        typedef mpl::true_ is_view;
        typedef mpl::int_<1> size;
        typedef T value_type;

        single_view()
            : val() {}

        explicit single_view(typename detail::call_param<T>::type val)
            : val(val) {}

        value_type val;
    };
    
    template <typename T>
    inline single_view<typename detail::as_fusion_element<T>::type>
    make_single_view(T const& v)
    {
        return single_view<typename detail::as_fusion_element<T>::type>(v);
    }
}}

#endif


