/*=============================================================================
    Copyright (c) 2001-2006 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_JOINT_VIEW_07162005_0140)
#define FUSION_JOINT_VIEW_07162005_0140

#include <boost/fusion/support/detail/access.hpp>
#include <boost/fusion/support/is_view.hpp>
#include <boost/fusion/sequence/intrinsic/begin.hpp>
#include <boost/fusion/sequence/intrinsic/end.hpp>
#include <boost/fusion/sequence/intrinsic/size.hpp>
#include <boost/fusion/sequence/view/joint_view/joint_view_iterator.hpp>
#include <boost/fusion/sequence/view/joint_view/detail/begin_impl.hpp>
#include <boost/fusion/sequence/view/joint_view/detail/end_impl.hpp>
#include <boost/fusion/support/sequence_base.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/bool.hpp>

namespace boost { namespace fusion
{
    struct joint_view_tag;
    struct forward_sequence_tag;
    struct fusion_sequence_tag;

    template <typename Sequence1, typename Sequence2>
    struct joint_view : sequence_base<joint_view<Sequence1, Sequence2> >
    {
        typedef joint_view_tag ftag;
        typedef fusion_sequence_tag tag; // this gets picked up by MPL
        typedef forward_sequence_tag category;
        typedef mpl::true_ is_view;

        typedef typename result_of::begin<Sequence1>::type first_type;
        typedef typename result_of::end<Sequence1>::type last_type;
        typedef typename result_of::begin<Sequence2>::type concat_type;
        typedef typename result_of::end<Sequence2>::type concat_last_type;
        typedef typename mpl::plus<result_of::size<Sequence1>, result_of::size<Sequence2> >::type size;

        joint_view(Sequence1& seq1, Sequence2& seq2)
            : seq1(seq1)
            , seq2(seq2)
        {}

        first_type first() const { return fusion::begin(seq1); }
        concat_type concat() const { return fusion::begin(seq2); }
        concat_last_type concat_last() const { return fusion::end(seq2); }

    private:

        typename mpl::if_<traits::is_view<Sequence1>, Sequence1, Sequence1&>::type seq1;
        typename mpl::if_<traits::is_view<Sequence2>, Sequence2, Sequence2&>::type seq2;
    };
}}

#endif


