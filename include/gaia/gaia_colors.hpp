//===================================================================
//
// Authors: Igor Chollet, Xavier Claeys, Pierre Fortin, Laura Grigori
//
//  Copyright (C) 2020 Sorbonne Universite, Inria
//  All right reserved.
//
//  This file is part of defmm.
//
//  defmm is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  defmm is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//  (see ../../LICENSE.txt)
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with defmm.  If not, see <http://www.gnu.org/licenses/>
//
//====================================================================

#ifndef GAIA_COLORS_HPP
#define GAIA_COLORS_HPP

#include <ostream>

namespace gaia{
  namespace colors{
      enum Code {
	FG_RED        = 31,
	FG_GREEN      = 32,
	FG_BLUE       = 34,
	FG_DEFAULT    = 39,
	FG_BLACK      = 30,
	FG_YELLOW     = 33,
	FG_MAGENTA    = 35,
	FG_CYAN       = 36,
	FG_LIGHTGREY  = 37,
	FG_DARKGREY   = 90,
	FG_WHITE      = 97,
	BG_RED        = 41,
	BG_GREEN      = 42,
	BG_BLUE       = 44,
	BG_DEFAULT    = 49,
	TXT_BOLD      = 1,
	TXT_DIM       = 2,
	TXT_UNDERLINE = 4,
	TXT_BLINK     = 5,
	TXT_DEFAULT   = 0
      };
    class metamorphe {
      Code code;
    public:
    metamorphe(Code pCode) : code(pCode) {}
      friend std::ostream&
      operator<<(std::ostream& os, const metamorphe& mod) {
	return os << "\033[" << mod.code << "m";
      }
    };
  }
} // GAIA

gaia::colors::metamorphe f_red      (gaia::colors::FG_RED       );
gaia::colors::metamorphe f_blue     (gaia::colors::FG_BLUE      );
gaia::colors::metamorphe f_green    (gaia::colors::FG_GREEN     );
gaia::colors::metamorphe f_def      (gaia::colors::FG_DEFAULT   );
gaia::colors::metamorphe f_black    (gaia::colors::FG_BLACK     );
gaia::colors::metamorphe f_yellow   (gaia::colors::FG_YELLOW    );
gaia::colors::metamorphe f_magenta  (gaia::colors::FG_MAGENTA   );
gaia::colors::metamorphe f_cyan     (gaia::colors::FG_CYAN      );
gaia::colors::metamorphe f_lgrey    (gaia::colors::FG_LIGHTGREY );
gaia::colors::metamorphe f_dgrey    (gaia::colors::FG_DARKGREY  );
gaia::colors::metamorphe f_white    (gaia::colors::FG_WHITE     );
gaia::colors::metamorphe b_red      (gaia::colors::BG_RED       );
gaia::colors::metamorphe b_blue     (gaia::colors::BG_BLUE      );
gaia::colors::metamorphe b_green    (gaia::colors::BG_GREEN     );
gaia::colors::metamorphe b_def      (gaia::colors::BG_DEFAULT   );
gaia::colors::metamorphe t_bold     (gaia::colors::TXT_BOLD     );
gaia::colors::metamorphe t_dim      (gaia::colors::TXT_DIM      );
gaia::colors::metamorphe t_underline(gaia::colors::TXT_UNDERLINE);
gaia::colors::metamorphe t_blink    (gaia::colors::TXT_BLINK    );
gaia::colors::metamorphe t_def      (gaia::colors::TXT_DEFAULT  );

#endif
