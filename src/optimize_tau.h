// /////////////////////////////////////////////////////////////////////////////
//
// Constriant Newton's method for \tau (Collection-level topic mixtures)
//
// This file is part of clda.
//
// Copyright (c) 2016  Clint P. George
//
// clda is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// clda is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.
//
// /////////////////////////////////////////////////////////////////////////////


# include "utils.h"


extern
  mat optimize_tau(
      mat vi_tau,
      mat vi_tau_rho_ss,
      vec collection_size,
      unsigned int num_collections,
      unsigned int num_topics,
      double alpha_h,
      double gamma_h,
      unsigned int max_opt_iter = 10,
      double step_size = .5,
      unsigned int verbose = 2
  );
