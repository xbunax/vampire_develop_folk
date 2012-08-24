//====================================================================================================
//
//       				                    	Fields
//
//  			 		Subroutines to calculate fields for the hamiltonian
//	 
//									Version 1.0 R Evans 20/10/2008
//
//==================================================================================================== 
#include "atoms.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "demag.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vmpi.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

//========================
//function prototypes
//========================

int calculate_exchange_fields(const int,const int);
int calculate_anisotropy_fields(const int,const int);
int calculate_cubic_anisotropy_fields(const int,const int);
int calculate_applied_fields(const int,const int);
int calculate_thermal_fields(const int,const int);
int calculate_dipolar_fields(const int,const int);
void calculate_hamr_fields(const int,const int);
void calculate_fmr_fields(const int,const int);
void calculate_surface_anisotropy_fields(const int,const int);

int calculate_spin_fields(const int start_index,const int end_index){
	//======================================================
	// 		Subroutine to calculate spin dependent fields
	//
	//			Version 1.0 R Evans 20/10/2008
	//======================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "calculate_spin_fields has been called" << std::endl;}
	
	// Initialise Total Spin Fields to zero
	fill (atoms::x_total_spin_field_array.begin()+start_index,atoms::x_total_spin_field_array.begin()+end_index,0.0);
	fill (atoms::y_total_spin_field_array.begin()+start_index,atoms::y_total_spin_field_array.begin()+end_index,0.0);
	fill (atoms::z_total_spin_field_array.begin()+start_index,atoms::z_total_spin_field_array.begin()+end_index,0.0);

	// Exchange Fields
	if(sim::hamiltonian_simulation_flags[0]==1) calculate_exchange_fields(start_index,end_index);
	
	// Anisotropy Fields
	if(sim::UniaxialScalarAnisotropy || sim::TensorAnisotropy) calculate_anisotropy_fields(start_index,end_index);
	if(sim::CubicScalarAnisotropy) calculate_cubic_anisotropy_fields(start_index,end_index);
	//if(sim::hamiltonian_simulation_flags[1]==3) calculate_local_anis_fields();
	if(sim::surface_anisotropy==true) calculate_surface_anisotropy_fields(start_index,end_index);
	// Spin Dependent Extra Fields
	//if(sim::hamiltonian_simulation_flags[4]==1) calculate_??_fields();
	
	return 0;
}

int calculate_external_fields(const int start_index,const int end_index){
	//======================================================
	// 		Subroutine to calculate external fields
	//
	//			Version 1.0 R Evans 20/10/2008
	//======================================================
	//const int num_atoms = atoms::num_atoms;

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "calculate_external_fields has been called" << std::endl;}

	// Initialise Total External Fields to zero
	fill (atoms::x_total_external_field_array.begin()+start_index,atoms::x_total_external_field_array.begin()+end_index,0.0);
	fill (atoms::y_total_external_field_array.begin()+start_index,atoms::y_total_external_field_array.begin()+end_index,0.0);
	fill (atoms::z_total_external_field_array.begin()+start_index,atoms::z_total_external_field_array.begin()+end_index,0.0);
	
	if(sim::program==7) calculate_hamr_fields(start_index,end_index);
	else{
	
		// Thermal Fields
		if(sim::hamiltonian_simulation_flags[3]==1) calculate_thermal_fields(start_index,end_index);

		// Applied Fields
		if(sim::hamiltonian_simulation_flags[2]==1) calculate_applied_fields(start_index,end_index);

	}
	
	// FMR Fields
	if(sim::hamiltonian_simulation_flags[5]==1) calculate_fmr_fields(start_index,end_index);

	// Dipolar Fields
	if(sim::hamiltonian_simulation_flags[4]==1) calculate_dipolar_fields(start_index,end_index);
	
	return 0;
}

int calculate_exchange_fields(const int start_index,const int end_index){
	//======================================================
	// 		Subroutine to calculate exchange fields
	//
	//			Version 2.0 Richard Evans 08/09/2011
	//======================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "calculate_exchange_fields has been called" << std::endl;}

	// Use appropriate function for exchange calculation
	switch(atoms::exchange_type){
		case 0: // isotropic
			for(int atom=start_index;atom<end_index;atom++){
				register double Hx=0.0;
				register double Hy=0.0;
				register double Hz=0.0;
				const int start=atoms::neighbour_list_start_index[atom];
				const int end=atoms::neighbour_list_end_index[atom]+1;
				for(int nn=start;nn<end;nn++){
					const int natom = atoms::neighbour_list_array[nn];
					const double Jij=atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij;
					Hx -= Jij*atoms::x_spin_array[natom];
					Hy -= Jij*atoms::y_spin_array[natom];
					Hz -= Jij*atoms::z_spin_array[natom];
				}
				atoms::x_total_spin_field_array[atom] += Hx;
				atoms::y_total_spin_field_array[atom] += Hy;
				atoms::z_total_spin_field_array[atom] += Hz;
			}
			break;
		case 1: // vector
			for(int atom=start_index;atom<end_index;atom++){
				register double Hx=0.0;
				register double Hy=0.0;
				register double Hz=0.0;
				const int start=atoms::neighbour_list_start_index[atom];
				const int end=atoms::neighbour_list_end_index[atom]+1;
				for(int nn=start;nn<end;nn++){
					const int natom = atoms::neighbour_list_array[nn];
					const int iid = atoms::neighbour_interaction_type_array[nn]; // interaction id
					const double Jij[3]={atoms::v_exchange_list[iid].Jij[0],
												atoms::v_exchange_list[iid].Jij[1],
												atoms::v_exchange_list[iid].Jij[2]};
					
					Hx -= Jij[0]*atoms::x_spin_array[natom];
					Hy -= Jij[1]*atoms::y_spin_array[natom];
					Hz -= Jij[2]*atoms::z_spin_array[natom];
				}
				atoms::x_total_spin_field_array[atom] += Hx;
				atoms::y_total_spin_field_array[atom] += Hy;
				atoms::z_total_spin_field_array[atom] += Hz;
			}
			break;
		case 2: // tensor
			for(int atom=start_index;atom<end_index;atom++){
				register double Hx=0.0;
				register double Hy=0.0;
				register double Hz=0.0;
				const int start=atoms::neighbour_list_start_index[atom];
				const int end=atoms::neighbour_list_end_index[atom]+1;
				for(int nn=start;nn<end;nn++){
					const int natom = atoms::neighbour_list_array[nn];
					const int iid = atoms::neighbour_interaction_type_array[nn]; // interaction id
					const double Jij[3][3]={atoms::t_exchange_list[iid].Jij[0][0],
													atoms::t_exchange_list[iid].Jij[0][1],
													atoms::t_exchange_list[iid].Jij[0][2],

													atoms::t_exchange_list[iid].Jij[1][0],
													atoms::t_exchange_list[iid].Jij[1][1],
													atoms::t_exchange_list[iid].Jij[1][2],

													atoms::t_exchange_list[iid].Jij[2][0],
													atoms::t_exchange_list[iid].Jij[2][1],
													atoms::t_exchange_list[iid].Jij[2][2]};
					
					const double S[3]={atoms::x_spin_array[natom],atoms::y_spin_array[natom],atoms::z_spin_array[natom]};
					
					Hx -= (Jij[0][0]*S[0] + Jij[0][1]*S[1] +Jij[0][2]*S[2]);
					Hy -= (Jij[1][0]*S[0] + Jij[1][1]*S[1] +Jij[1][2]*S[2]);
					Hz -= (Jij[2][0]*S[0] + Jij[2][1]*S[1] +Jij[2][2]*S[2]);
				}
				atoms::x_total_spin_field_array[atom] += Hx;
				atoms::y_total_spin_field_array[atom] += Hy;
				atoms::z_total_spin_field_array[atom] += Hz;
			}
			break;
		}

		return EXIT_SUCCESS;
	}

int calculate_anisotropy_fields(const int start_index,const int end_index){
	//======================================================
	// 	Subroutine to calculate uniaxial anisotropy fields
	//
	//			Version 1.0 R Evans 20/10/2008
	//======================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "calculate_anisotropy_fields has been called" << std::endl;}

		// Use appropriate function for anisotropy calculation
	switch(sim::AnisotropyType){
		case 0: // scalar
			for(int atom=start_index;atom<end_index;atom++){
				const int imaterial=atoms::type_array[atom];
				atoms::z_total_spin_field_array[atom] -= 2.0*mp::MaterialScalarAnisotropyArray[imaterial].K*atoms::z_spin_array[atom];
			}
			break;
		case 1: // tensor
			for(int atom=start_index;atom<end_index;atom++){
				const int imaterial=atoms::type_array[atom];

				const double K[3][3]={2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[0][0],
												2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[0][1],
												2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[0][2],

												2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[1][0],
												2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[1][1],
												2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[1][2],

												2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[2][0],
												2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[2][1],
												2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[2][2]};
					
				const double S[3]={atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};

				atoms::x_total_spin_field_array[atom] -= (K[0][0]*S[0] + K[0][1]*S[1] +K[0][2]*S[2]);
				atoms::y_total_spin_field_array[atom] -= (K[1][0]*S[0] + K[1][1]*S[1] +K[1][2]*S[2]);
				atoms::z_total_spin_field_array[atom] -= (K[2][0]*S[0] + K[2][1]*S[1] +K[2][2]*S[2]);
			}
			break;
		}

	return EXIT_SUCCESS;
}

int calculate_cubic_anisotropy_fields(const int start_index,const int end_index){
	//------------------------------------------------------
	// 	Function to calculate cubic anisotropy fields
	//
	//			Version 1.0 R Evans 28/07/2012
	//
	//		E = -0.5 Kc (Sx^4 + Sy^4 + Sz^4)
	//		Hx = +2 Kc*(Sx^3)
	//		Hy = +2 Kc*(Sy^3)
	//		Hz = +2 Kc*(Sz^3)
	//	
	//------------------------------------------------------
	//std::cout << "here" << std::endl;
	for(int atom=start_index;atom<end_index;atom++){
		const int imaterial=atoms::type_array[atom];
		const double Kc=2.0*mp::MaterialCubicAnisotropyArray[imaterial];

		const double Sx=atoms::x_spin_array[atom];
		atoms::x_total_spin_field_array[atom] -= Kc*Sx*Sx*Sx;
		
		const double Sy=atoms::y_spin_array[atom];
		atoms::y_total_spin_field_array[atom] -= Kc*Sy*Sy*Sy;

		const double Sz=atoms::z_spin_array[atom];
		atoms::z_total_spin_field_array[atom] -= Kc*Sz*Sz*Sz;
		
	}
	return EXIT_SUCCESS;
}

void calculate_surface_anisotropy_fields(const int start_index,const int end_index){
	//======================================================
	// 		Subroutine to calculate surface anisotropy fields
	//
	//			Version 1.0 Richard Evans 13/09/2011
	//======================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "calculate_surface_anisotropy_fields has been called" << std::endl;}

	for(int atom=start_index;atom<end_index;atom++){
		// only calculate for surface atoms
		if(atoms::surface_array[atom]==true){
			const int imaterial=atoms::type_array[atom];
			const double Ks=2.0*mp::material[imaterial].Ks; // note factor two here from differentiation
			const double S[3]={atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
		
			for(int nn=atoms::nearest_neighbour_list_si[atom];nn<atoms::nearest_neighbour_list_ei[atom];nn++){
				const double si_dot_eij=(S[0]*atoms::eijx[nn]+S[1]*atoms::eijy[nn]+S[2]*atoms::eijz[nn]);
				atoms::x_total_spin_field_array[atom]-=Ks*si_dot_eij*atoms::eijx[nn];
				atoms::y_total_spin_field_array[atom]-=Ks*si_dot_eij*atoms::eijy[nn];
				atoms::z_total_spin_field_array[atom]-=Ks*si_dot_eij*atoms::eijz[nn];
			}
		}
	}
	
	return;
}

int calculate_applied_fields(const int start_index,const int end_index){
	//======================================================
	// 	Subroutine to calculate applied fields
	//
	//			Version 1.0 R Evans 20/10/2008
	//======================================================
	//const int num_atoms = atoms::num_atoms;

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "calculate_applied_fields has been called" << std::endl;}

	for(int atom=start_index;atom<end_index;atom++){
		atoms::x_total_external_field_array[atom] += sim::H_vec[0]*sim::H_applied;
		atoms::y_total_external_field_array[atom] += sim::H_vec[1]*sim::H_applied;
		atoms::z_total_external_field_array[atom] += sim::H_vec[2]*sim::H_applied;

		//std::cout << atom << "\tapplied fields\t" << sim::H_vec[0]*sim::H_applied << "\t";
		//std::cout << sim::H_vec[1]*sim::H_applied << "\t";
		//std::cout << sim::H_vec[2]*sim::H_applied << std::endl;
	}

	// Add external field from thin film sample
	if(sim::ext_demag==true){
		
		// calculate system magnetisation
		stats::mag_m();
		
		// calculate global demag field -mu_0 M D, M = m/V
		const double mu_0= -4.0*M_PI*1.0e-7/(cs::system_dimensions[0]*cs::system_dimensions[1]*cs::system_dimensions[2]*1.0e-30);
		const double HD[3]={	mu_0*sim::demag_factor[0]*stats::total_mag_actual[0],
									mu_0*sim::demag_factor[1]*stats::total_mag_actual[1],
									mu_0*sim::demag_factor[2]*stats::total_mag_actual[2]};
		
		//std::cout << "mu_0" << "\t" << mu_0 << std::endl;
		//std::cout << "Magnetisation " << stats::total_mag_actual[0] << "\t" << stats::total_mag_actual[1] << "\t" << stats::total_mag_actual[2] << std::endl;  
		//std::cout << "External Demag Field " << HD[0] << "\t" << HD[1] << "\t" << HD[2] << std::endl;  
		for(int atom=start_index;atom<end_index;atom++){
			atoms::x_total_external_field_array[atom] += HD[0];
			atoms::y_total_external_field_array[atom] += HD[1];
			atoms::z_total_external_field_array[atom] += HD[2];
		}
	}

	return 0;
}

int calculate_thermal_fields(const int start_index,const int end_index){
	//======================================================
	// 		Subroutine to calculate thermal fields
	//
	//			Version 1.1 R Evans 26/07/2012
	//======================================================

	const double sqrt_T=sqrt(sim::temperature);

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "calculate_thermal_fields has been called" << std::endl;}

	// unroll Sigma
	std::vector<double> SigmaPre(0);
	for(int mat=0;mat<mp::material.size();mat++) SigmaPre.push_back(sqrt_T*mp::material[mat].H_th_sigma);

 	generate (atoms::x_total_external_field_array.begin()+start_index,atoms::x_total_external_field_array.begin()+end_index, mtrandom::gaussian);
	generate (atoms::y_total_external_field_array.begin()+start_index,atoms::y_total_external_field_array.begin()+end_index, mtrandom::gaussian);
	generate (atoms::z_total_external_field_array.begin()+start_index,atoms::z_total_external_field_array.begin()+end_index, mtrandom::gaussian);

	for(int atom=start_index;atom<end_index;atom++){

		const int imaterial=atoms::type_array[atom];
		const double H_th_sigma = SigmaPre[imaterial];

		atoms::x_total_external_field_array[atom] *= H_th_sigma;
		atoms::y_total_external_field_array[atom] *= H_th_sigma;
		atoms::z_total_external_field_array[atom] *= H_th_sigma; 
	}

	return EXIT_SUCCESS;
}

int calculate_dipolar_fields(const int start_index,const int end_index){
	//======================================================
	// 		Subroutine to calculate dipolar fields
	//
	//			Version 1.0 R Evans 02/11/2009
	//======================================================

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "calculate_dipolar_fields has been called" << std::endl;}

	// Add dipolar fields
	for(int atom=start_index;atom<end_index;atom++){
		atoms::x_total_external_field_array[atom] += atoms::x_dipolar_field_array[atom];
		atoms::y_total_external_field_array[atom] += atoms::y_dipolar_field_array[atom];
		atoms::z_total_external_field_array[atom] += atoms::z_dipolar_field_array[atom];
	}

	return 0;
}

void calculate_hamr_fields(const int start_index,const int end_index){
	
	if(err::check==true){std::cout << "calculate_hamr_fields has been called" << std::endl;}

	// Declare hamr variables
	const double fwhm=200.0; // A
	const double fwhm2=fwhm*fwhm;
	const double px = sim::head_position[0];
	const double py = sim::head_position[1];
	const double DeltaT=sim::Tmax-sim::Tmin;

	// declare head-field variables
	const double H_bounds_min[2]={-400.0,-250.0}; // A
	const double H_bounds_max[2]={-100.0,+250.0}; // A
	const double H_osc_freq=200.0; // A
	const double Hloc_min_x=sim::head_position[0]+H_bounds_min[0];
	const double Hloc_min_y=sim::head_position[1]+H_bounds_min[1];
	const double Hloc_max_x=sim::head_position[0]+H_bounds_max[0];
	const double Hloc_max_y=sim::head_position[1]+H_bounds_max[1];
	const double Hloc_parity_field=sim::H_applied*double(2*(int(sim::head_position[0]/H_osc_freq)%2)-1);
	const double Hvecx=sim::H_vec[0];
	const double Hvecy=sim::H_vec[1];
	const double Hvecz=sim::H_vec[2];
	
	// Add localised thermal field
	generate (atoms::x_total_external_field_array.begin()+start_index,atoms::x_total_external_field_array.begin()+end_index, mtrandom::gaussian);
	generate (atoms::y_total_external_field_array.begin()+start_index,atoms::y_total_external_field_array.begin()+end_index, mtrandom::gaussian);
	generate (atoms::z_total_external_field_array.begin()+start_index,atoms::z_total_external_field_array.begin()+end_index, mtrandom::gaussian);

	if(sim::head_laser_on){
		for(int atom=start_index;atom<end_index;atom++){
			const int imaterial=atoms::type_array[atom];
			const double cx = atoms::x_coord_array[atom];
			const double cy = atoms::y_coord_array[atom];		
			const double r2 = (cx-px)*(cx-px)+(cy-py)*(cy-py);
			const double sqrt_T = sqrt(sim::Tmin+DeltaT*exp(-r2/fwhm2));
			const double H_th_sigma = sqrt_T*mp::material[imaterial].H_th_sigma;
			atoms::x_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
			atoms::y_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
			atoms::z_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
		}

		// Add localised applied field
		for(int atom=start_index;atom<end_index;atom++){
			const double cx = atoms::x_coord_array[atom];
			const double cy = atoms::y_coord_array[atom];		
			double Hx=0.0;
			double Hy=0.0;
			double Hz=0.0;
			if((cx >= Hloc_min_x) && (cx <= Hloc_max_x) && (cy >= Hloc_min_y) && (cy <= Hloc_max_y)){
				Hx=Hvecx*Hloc_parity_field;
				Hy=Hvecy*Hloc_parity_field;
				Hz=Hvecz*Hloc_parity_field;
			}
			atoms::x_total_external_field_array[atom] += Hx;
			atoms::y_total_external_field_array[atom] += Hy;
			atoms::z_total_external_field_array[atom] += Hz;
		}
	}
	else{
		// Otherwise just use global temperature
		double sqrt_T=sqrt(sim::temperature);
		for(int atom=start_index;atom<end_index;atom++){
			const int imaterial=atoms::type_array[atom];
			const double H_th_sigma = sqrt_T*material_parameters::material[imaterial].H_th_sigma;
			atoms::x_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
			atoms::y_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
			atoms::z_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
		}
	}
}

void calculate_fmr_fields(const int start_index,const int end_index){
	
	if(err::check==true){std::cout << "calculate_fmr_fields has been called" << std::endl;}

	// Declare fmr variables
	const double real_time=sim::time*mp::dt_SI;
	const double osc_freq=20.0e9; // Hz
	const double osc_period=1.0/osc_freq;
	const double Hfmr_vec[3]={1.0,0.0,0.0};
	const double Hfmr=0.001; // T
	const double Hx=Hfmr_vec[0]*Hfmr*sin(2.0*M_PI*real_time/osc_period);
	const double Hy=Hfmr_vec[1]*Hfmr*sin(2.0*M_PI*real_time/osc_period);
	const double Hz=Hfmr_vec[2]*Hfmr*sin(2.0*M_PI*real_time/osc_period);
	
	// Add localised applied field
	for(int atom=start_index;atom<end_index;atom++){
			atoms::x_total_external_field_array[atom] += Hx;
			atoms::y_total_external_field_array[atom] += Hy;
			atoms::z_total_external_field_array[atom] += Hz;
	}
}

