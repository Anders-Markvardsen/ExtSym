#include "data_types.h"  // description declared in data_types.h
#include <string.h>
#include <stdlib.h>  // for exit(1)
#include <stdio.h>

struct description *symmetries;

int num_extinction_groups;		// the number of extinction groups


struct description ortho_symmetries[] ={
        {"", ""},     
	{"________",	"P -  -  -  "}, 
	{"______l_",	"P -  -  21 "},
	{"_____k__",	"P -  21 -  "},
	{"_____kl_",	"P -  21 21 "},
	{"____h___",	"P 21 -  -  "},
	{"____h_l_",	"P 21 -  21 "},
	{"____hk__",	"P 21 21 -  "},
	{"____hkl_",	"P 21 21 21 "},
	{"___h____",	"P -  -  a  "},  // hk0 with h=odd absent implies h00 with h=odd absent
	{"___k____",	"P -  -  b  "},   
	{"___s____",	"P -  -  n  "},  
	{"__h_____",	"P -  a  -  "},
	{"__hh____",	"P -  a  a  "},
	{"__hk____",	"P -  a  b  "},
	{"__hs____",	"P -  a  n  "},
	{"__l_____",	"P -  c  -  "},
	{"__lh____",	"P -  c  a  "},
	{"__lk____",	"P -  c  b  "},
	{"__ls____",	"P -  c  n  "}, 
	{"__s_____",	"P -  n  -  "},// 20 extinction group
	{"__sh____",	"P -  n  a  "},
	{"__sk____",	"P -  n  b  "},
	{"__ss____",	"P -  n  n  "},
	{"_k______",	"P b  -  -  "},
	{"_k_h____",	"P b  -  a  "},
	{"_k_k____",	"P b  -  b  "},
	{"_k_s____",	"P b  -  n  "},
	{"_kh_____",	"P b  a  -  "},
	{"_khh____",	"P b  a  a  "},
	{"_khk____",	"P b  a  b  "},
	{"_khs____",	"P b  a  n  "},
	{"_kl_____",	"P b  c  -  "},
	{"_klh____",	"P b  c  a  "},
	{"_klk____",	"P b  c  b  "},
	{"_kls____",	"P b  c  n  "},
	{"_ks_____",	"P b  n  -  "},
	{"_ksh____",	"P b  n  a  "},
	{"_ksk____",	"P b  n  b  "},
	{"_kss____",	"P b  n  n  "},
	{"_l______",	"P c  -  -  "}, // 40
	{"_l_h____",	"P c  -  a  "},
	{"_l_k____",	"P c  -  b  "},
	{"_l_s____",	"P c  -  n  "},
	{"_lh_____",	"P c  a  -  "},
	{"_lhh____",	"P c  a  a  "},
	{"_lhk____",	"P c  a  b  "},
	{"_lhs____",	"P c  a  n  "},
	{"_ll_____",	"P c  c  -  "},
	{"_llh____",	"P c  c  a  "},
	{"_llk____",	"P c  c  b  "},
	{"_lls____",	"P c  c  n  "},
	{"_ls_____",	"P c  n  -  "},
	{"_lsh____",	"P c  n  a  "},
	{"_lsk____",	"P c  n  b  "},
	{"_lss____",	"P c  n  n  "},
	{"_s______",	"P n  -  -  "},
	{"_s_h____",	"P n  -  a  "},
	{"_s_k____",	"P n  -  b  "},
	{"_s_s____",	"P n  -  n  "},
	{"_sh_____",	"P n  a  -  "}, // 60
	{"_shh____",	"P n  a  a  "},
	{"_shk____",	"P n  a  b  "},
	{"_shs____",	"P n  a  n  "},
	{"_sl_____",	"P n  c  -  "},
	{"_slh____",	"P n  c  a  "},
	{"_slk____",	"P n  c  b  "},
	{"_sls____",	"P n  c  n  "},
	{"_ss_____",	"P n  n  -  "},
	{"_ssh____",	"P n  n  a  "},
	{"_ssk____",	"P n  n  b  "},
	{"_sss____",	"P n  n  n  "},   // corrected from _ssk___ to _sss___
	{"C_______",	"C -  -  -  "},
	{"C_____l_",	"C -  -  21 "},
	{"C__h____",	"C -  - (ab)"},
	{"C_l_____",	"C -  c  -  "},
	{"C_lh____",	"C -  c (ab)"},
	{"Cl______",	"C c  -  -  "},
	{"Cl_h____",	"C c  - (ab)"},
	{"Cll_____",	"C c  c  -  "},
	{"Cllh____",	"C c  c (ab)"}, // 80
	{"B_______",	"B -  -  -  "},
	{"B____k__",	"B -  21 -  "},
	{"B__k____",	"B -  -  b  "},
	{"B_h_____",	"B - (ac)-  "},
	{"B_hk____",	"B - (ac)b  "},
	{"Bk______",	"B b  -  -  "},
	{"Bk_k____",	"B b  -  b  "},
	{"Bkh_____",	"B b (ac)-  "},
	{"Bkhk____",	"B b (ac)b  "},
	{"A_______",	"A -  -  -  "},
	{"A___h___",	"A 21 -  -  "},
	{"A__h____",	"A -  -  a  "},
	{"A_h_____",	"A -  a  -  "},
	{"A_hh____",	"A -  a  a  "},
	{"Ak______",	"A(bc)-  -  "},
	{"Ak_h____",	"A(bc)-  a  "},
	{"Akh_____",	"A(bc)a  -  "},
	{"Akhh____",	"A(bc)a  a  "},
	{"I_______",	"I -  -  -  "},
	{"I__h____",	"I -  - (ab)"}, // 100
	{"I_h_____",	"I - (ac)-  "},
	{"I_hh____",	"I -  c  b  "},
	{"Ik______",	"I(bc)-  -  "},
	{"Ik_h____",	"I c  -  a  "},
	{"Ikh_____",	"I b  a  -  "},
	{"Ikhh____",	"I b  c  a  "},

	{"F_______", "F -  -  -  "},

	{"Fk44____", "F -  d  d  "},	// The 4-symbol refer to: in addition to miller

	{"F4h4____", "F d  -  d  "}, // indices being odd there sum must also be also be

	{"F44h____", "F d  d  -  "}, // different from a multiple of 4 to be a present

	{"F444____", "F d  d  d  "}  // intensity.
};    


struct description monob_symmetries[] ={
        {"", ""},  
	{"________",	"P 1  -  1  "}, 
	{"_____k__",	"P 1  21 1  "},
	{"__h_____",	"P 1  a  1  "},
	{"__h__k__",	"P 1 21/a 1 "},
	{"__l_____",	"P 1  c  1  "},
	{"__l__k__",	"P 1 21/c 1 "},
	{"__s_____",	"P 1  n  1  "},
	{"__s__k__",	"P 1 21/n 1 "},
	{"C_______",	"C 1  -  1  "},
	{"C_l_____",	"C 1  c  1  "},
	{"A_______",	"A 1  -  1  "},
	{"A_h_____",	"A 1  n  1  "},
	{"I_______",	"I 1  -  1  "},
	{"I_h_____",	"I 1  a  1  "},
};

struct description monoc_symmetries[] ={
        {"", ""},  
	{"________",	"P 1  1  -  "}, 
	{"______l_",	"P 1  1  21 "},
	{"___h____",	"P 1  1  a  "},
	{"___h__l_",	"P 1  1 21/a"},
	{"___k____",	"P 1  1  b  "},
	{"___k__l_",	"P 1  1 21/c"},
	{"___s____",	"P 1  1  n  "},
	{"___s__l_",	"P 1  1 21/n"},
	{"B_______",	"B 1  1  -  "},
	{"B__k____",	"B 1  1  n  "},
	{"A_______",	"A 1  1  -  "},
	{"A__h____",	"A 1  1  a  "},
	{"I_______",	"I 1  1  -  "},
	{"I__h____",	"I 1  1  b  "},
};

struct description monoa_symmetries[] ={
        {"", ""},  
	{"________",	"P -  1  1  "}, 
	{"____h___",	"P 21 1  1  "},
	{"_k______",	"P b  1  1  "},
	{"_k__h___",	"P 21/b 1 1 "},
	{"_l______",	"P c  1  1  "},
	{"_l__h___",	"P 21/c 1 1 "},
	{"_s______",	"P n  1  1  "},
	{"_s__h___",	"P 21/n 1 1 "},
	{"C_______",	"C -  1  1  "},
	{"Cl______",	"C n  1  1  "},
	{"B_______",	"B -  1  1  "},
	{"Bk______",	"B b  1  1  "},
	{"I_______",	"I -  1  1  "},
	{"Ik______",	"I c  1  1  "}
};

struct description tetragonal_symmetries[] ={
        {"", ""},  
	{"________",	"P -  -  -  "}, 
	{"____hk__",	"P -  21 -  "},
	{"______l_",	"P 42 -  -  "},
	{"____hkl_",	"P 42 21 -  "},
	{"______4_",	"P 41 -  -  "},
	{"____hk4_",	"P 41 21 -  "},
	{"_______l",	"P -  -  c  "},
	{"____hk_l",	"P -  21 c  "},
	{"_kh_____",	"P -  b  -  "},
	{"_kh____l",	"P -  b  c  "},
	{"_ll_____",	"P -  c  -  "},
	{"_ll____l",	"P -  c  c  "},
	{"_ss_____",	"P -  n  -  "},
	{"_ss____l",	"P -  n  c  "},
	{"___s____",	"P n  -  -  "}, 
	{"___s__l_",	"P 42/n - - "},
	{"___s___l",	"P n  -  c  "},
	{"_khs____",	"P n  b  -  "},
	{"_khs___l",	"P n  b  c  "},
	{"_lls____",	"P n  c  -  "},
	{"_lls___l",	"P n  c  c  "},
	{"_sss____",	"P n  n  -  "},
	{"_sss___l",	"P n  n  c  "},
	{"I_______",	"I -  -  -  "},
	{"I_____4_",	"I 41 -  -  "},
	{"I______4",	"I -  -  d  "},
	{"Ikh_____",	"I -  c  -  "},
	{"Ikh____4",	"I -  c  d  "},
	{"I__h__4_",	"I 41/a - - "},
	{"I__h___4",	"I a  -  d  "},
	{"Ikhh___4",	"I a  c  d  "}
};


struct description trigonal_symmetries[] ={  
        {"", ""},  
	{"________",	"P -  -  -  "}, 
	{"______3_",	"P 31  - -  "},
	{"_______l",	"P -  -  c  "},  //(h=h)
	{"_______i",	"P -  c  -  "},  //(h=-h)
	{"R_______",	"R -  -  -  "},  
	{"R______l",	"R -  -  c  "}
};


struct description hexagonal_symmetries[] ={  
        {"", ""},  
	{"________",	"P -  -  -  "}, 
	{"______l_",	"P 63 -  -  "},
	{"______3_",	"P 62 -  -  "},  
	{"______6_",	"P 61 -  -  "},  
	{"_______l",	"P -  -  c  "},  //(h=h)
	{"_______i",	"P -  c  -  "},  //(h=-h)
	{"_______e",	"P -  c  c  "}		// both
};


struct description tri_and_hexa_symmetries[] ={  // hexagonal axes and reverse setting
        {"", ""},  
	{"________",	"P -  -  -  "}, 
	{"______l_",	"P 63 -  -  "},
	{"______3_",	"P 31  - -  "},
	{"______3_",	"P 62 -  -  "}, 
	{"______6_",	"P 61 -  -  "},
	{"_______l",	"P -  -  c  "}, 
	{"_______i",	"P -  c  -  "},
	{"_______e",	"P -  c  c  "},	
	{"R_______",	"R -  -  -  "},  
	{"R______l",	"R -  -  c  "}
};

struct description cubic_symmetries[] ={
        {"", ""},  
	{"________",	"P -  -  -  "}, 
	{"____hkl_",	"P 21(42) --"},
	{"____444_",	"P 41 -  -  "},
	{"_______c",	"P -  -  n  "},
	{"_klh____",	"P a  -  -  "},
	{"_sss____",	"P n  -  -  "},
	{"_sss___c",	"P n  -  n  "},
	{"I_______",	"I -  -  -  "},
	{"I___444_",	"I 41 -  -  "},
	{"I______b",	"I -  -  d  "},
	{"Iklh____",	"I a  -  -  "},
	{"Iklh___b",	"I a  -  d  "},
	{"F_______",	"F -  -  -  "},
	{"F___444_",	"F 41 -  -  "},
	{"F______c",	"F -  -  c  "}, 
	{"F444____",	"F d  -  -  "},
	{"F444___c",	"F d  -  c  "},
};

void initialize_sym_struct(char laue_class_symbol[])
{
	if ( laue_class_symbol[0] == 'o' ) {
			num_extinction_groups = 111;
			symmetries = ortho_symmetries;
	}
	else if ( laue_class_symbol[0] == 'b') {
			num_extinction_groups = 14;
			symmetries = monob_symmetries;
	}
	else if ( laue_class_symbol[0] == 'c') {
			num_extinction_groups = 14;
			symmetries = monoc_symmetries;
	}
	else if ( laue_class_symbol[0] == 'a') {
			num_extinction_groups = 14;
			symmetries = monoa_symmetries;
	}
	else if ( laue_class_symbol[0] == 't') {
			num_extinction_groups = 31;
			symmetries = tetragonal_symmetries;
	}
	else if ( laue_class_symbol[0] == 'r') {
			num_extinction_groups = 6;
			symmetries = trigonal_symmetries;
	}
	else if ( laue_class_symbol[0] == 'h') {
			num_extinction_groups = 7;
			symmetries = hexagonal_symmetries;
	}
	else if ( laue_class_symbol[0] == 'd') {
			num_extinction_groups = 17;
			symmetries = cubic_symmetries;
	}
	else if ( laue_class_symbol[0] == 'x') {
			num_extinction_groups = 10;
			symmetries = tri_and_hexa_symmetries;
	}
	else {
			printf("laue_class_symbol not valid\n");
			exit(1);
	}
}





