import json
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
np.set_printoptions(precision=10)

class finite_wing: #### the code for W_1 is assuming a symmetric airfoil (alpha_L0 = 0.0)
    """This class contains functions that calculates position of nodes, control points, li, xi, eta, phi, psi, P_matrix, A_matrix, gamma_vector, cartesian_velocity, C_p, C_L, C_mle, C_mc/4"""
    def __init__(self, json_file):
        """This initializes the finite wing class"""
        self.json_file = json_file
        self.load_json()

    def load_json(self):
        """This function pulls in all the input values from the json"""
        with open(self.json_file, 'r') as json_handle:
            input_vals = json.load(json_handle)
            self.wing_type = input_vals["wing"]["planform"]["type"]
            self.aspect_ratio = input_vals["wing"]["planform"]["aspect_ratio"]
            self.taper_ratio = input_vals["wing"]["planform"]["taper_ratio"]
            self.airfoil_lift_slope = input_vals["wing"]["airfoil_lift_slope"]
            self.nodes_per_semispan = input_vals["wing"]["nodes_per_semispan"]
            self.washout_distribution = input_vals["wing"]["washout"]["distribution"]
            self.Omega = input_vals["wing"]["washout"]["amount[deg]"]
            self.CL_design = input_vals["wing"]["washout"]["CL_design"]
            self.begin_first_aileron_zb = input_vals["wing"]["aileron"]["begin[z/b]"]
            self.end_first_aileron_zb = input_vals["wing"]["aileron"]["end[z/b]"]
            self.begin_aileron_cfc = input_vals["wing"]["aileron"]["begin[cf/c]"]
            self.end_aileron_cfc = input_vals["wing"]["aileron"]["end[cf/c]"]
            self.begin_second_aileron_zb = -input_vals["wing"]["aileron"]["begin[z/b]"]
            self.end_second_aileron_zb = -input_vals["wing"]["aileron"]["end[z/b]"]
            self.hinge_efficiency = input_vals["wing"]["aileron"]["hinge_efficiency"]
            self.deflection_efficiency = input_vals["wing"]["aileron"]["deflection_efficiency"]
            self.N_nodes = (2* self.nodes_per_semispan) - 1
            # self.alphas = input_vals["condition"]["alpha_root[deg]"]
            self.alpha_root = np.radians(input_vals["condition"]["alpha_root[deg]"])
            self.aileron_deflection = np.radians(input_vals["condition"]["aileron_deflection[deg]"])
            self.pbar = input_vals["condition"]['pbar']
            self.planform = input_vals["view"]["planform"]
            self.is_washout_distribution = input_vals["view"]["washout_distribution"]
            self.is_aileron_distribution = input_vals["view"]["aileron_distribution"]


if __name__ == "__main__":
    wing = finite_wing("W3_input.json")
    print("\nalpha!!!!!!!!!!!!!!!!:", np.degrees(wing.alpha_root), "\n")
    # first calculate array of thetas (should be the size of N_nodes)
    theta_array = np.zeros((wing.N_nodes))
    # calculate the first theta value
    theta_array[0] = 0
    # calculate the last theta value
    theta_array[wing.N_nodes - 1] = np.pi
    #calculate the middle theta values
    for i in range(1, wing.N_nodes-1):
        theta_array[i] = (((i)*np.pi)/(wing.N_nodes-1))
    print("\n") 
    print("theta")  
    print(theta_array)
    print("\n")

    # calculate array of z/b    
    z_b_array = np.zeros((wing.N_nodes))
    z_b_array[0] = 0.5
    z_b_array[wing.N_nodes - 1] = -0.5
    # now, calculate the middle z values
    for i in range(1, wing.N_nodes -1):
        z_b_array[i] = (1/2)*np.cos(theta_array[i]) # 1/2 because b = 1
    print("z_b_array")
    print(z_b_array)
    print("\n")

    # calculate chord distribution array
    chord_dist_array = np.zeros((wing.N_nodes))
    for i in range(0, wing.N_nodes):
        if wing.wing_type == "elliptic":
            chord_dist_array[i] = ((4*1)/(np.pi*wing.aspect_ratio))*np.sin(theta_array[i])
        elif wing.wing_type == "tapered":
            chord_dist_array[i] = (2*1/(wing.aspect_ratio*(1 + wing.taper_ratio)))*(1-((1-wing.taper_ratio)*abs(np.cos(theta_array[i]))))
    print("chord_distribution")
    print(chord_dist_array)
    print("\n")

    # now we have the requisite information to fill up the C_matrix 
    C_matrix = np.zeros((wing.N_nodes, wing.N_nodes))

    # first, fill in the first row. 
    for j in range(0, wing.N_nodes + 1):
        C_matrix[0, j-1] = (j)**2

    # now, fill in the last row
    for j in range(0, wing.N_nodes+1):
        C_matrix[wing.N_nodes -1, j-1] = ((-1) ** (j + 1)) * (j ** 2)
        # C_matrix[, j-1] = (j)**2
    # print("C_matrix \n", C_matrix, "\n")

    for i in range(1, wing.N_nodes-1):
        for j in range(0, wing.N_nodes):
            C_matrix[i, j] = (((4 * 1) / (wing.airfoil_lift_slope * chord_dist_array[i])) + (j+1) / np.sin(theta_array[i])) * np.sin((j+1) * theta_array[i])


    # now, calculate the inverse C_matrix
    C_matrix_inv = np.linalg.inv(C_matrix)


    # create the B vector (should all be ones to find the fourier coefficients)
    B_vector = np.zeros((wing.N_nodes))
    for i in range(0, wing.N_nodes):
        B_vector[i] = 1.0
    # print("B_vector:\n", B_vector, "\n")

    # now multiply the C inverse matrix by the B vector to get the fourier coefficients (a_vector)
    a_vector = np.matmul(C_matrix_inv, B_vector)


    ####Now, get the kappa L value (only applies to tapered or square wings, but won't affect elliptic answer)
    if wing.wing_type == "elliptic":
        kappa_L = 0.0
    elif wing.wing_type == "tapered":
        kappa_L = (1-(1+(np.pi*wing.aspect_ratio/wing.airfoil_lift_slope))*a_vector[0])/((1+(np.pi*wing.aspect_ratio/wing.airfoil_lift_slope))*a_vector[0])
    print("\nKappa_L:\n", kappa_L, "\n")


    ####Now, calculate the kappa_D value (only applies to tapered or square wings, but won't affect elliptic answer)
    if wing.wing_type == "elliptic":
        kappa_D = 0.0
    elif wing.wing_type == "tapered":
        kappa_D = 0.0
        for i in range(1, len(a_vector)):
            kappa_D = kappa_D + (i+1)*((a_vector[i]/a_vector[0])**2)
    print("Kappa_D:\n", kappa_D, "\n")


    #Now, calculate the e_s value (only applies to tapered or square wings, but won't affect elliptic answer)
    e_s = 1/(1 + kappa_D)
    print("e_s:\n", e_s, "\n")


    #### now, calculate the washout distribution. 
    if wing.is_washout_distribution:
        washout_distribution = np.zeros((wing.N_nodes))
        if wing.washout_distribution == "optimum":
            for i in range(0, wing.N_nodes):
                if (i == 0 or i == wing.N_nodes) and (wing.wing_type == "elliptic"):
                    washout_distribution[i] = 0.0
                else:
                    washout_distribution[i] = 1 - ((np.sin(theta_array[i]))/((chord_dist_array[i])/chord_dist_array[len(chord_dist_array) // 2])) # eq 6.31 in eng handbook
        elif wing.washout_distribution == "linear":
            for i in range(0, wing.N_nodes):
                washout_distribution[i] = abs(np.cos(theta_array[i])) # eq 6.30 in eng handbook
        elif wing.washout_distribution == "none":
            for i in range(0, wing.N_nodes):
                washout_distribution[i] = 0.0
            
    # print("\nwashout_distribution")
    # print(washout_distribution)
    # print("\n")   


    #### now, calculate the bn array by multiplying the C_inverse matrix by the values from the washout distribution
    b_vector = np.matmul(C_matrix_inv, washout_distribution)
    # add the b vector to the solutions text file along with the other stuff from before


    #### now, calculate epsilon_omega
    epsilon_omega = b_vector[0]/a_vector[0] # eq 6.22 aero eng handbook
    print("epsilon_omega:")
    print(epsilon_omega, "\n")


    #### now, calculate kappa_DL (same equations for elliptic and tapered)
    kappa_DL = 0.0
    for i in range(1, len(a_vector)):
        kappa_DL = kappa_DL + (i+1)*(a_vector[i]/a_vector[0])*((b_vector[i]/b_vector[0])-(a_vector[i]/a_vector[0])) # equation 1.8.30 mech of flight
    kappa_DL = 2*(b_vector[0]/a_vector[0])*kappa_DL
    print("kappa_DL:")
    print(kappa_DL,"\n")


    #### now, calculate kappa_DOmega (same equations for elliptic and tapered)
    kappa_DOmega = 0.0
    for i in range(1, len(a_vector)):
        kappa_DOmega = kappa_DOmega + (i+1)*(((b_vector[i]/b_vector[0])-(a_vector[i]/a_vector[0]))**2) # equation 1.8.31 mech of flight
    kappa_DOmega = ((b_vector[0]/a_vector[0])**2)*kappa_DOmega
    print("kappa_DOmega:")
    print(kappa_DOmega,"\n")

    #### now, re-evaluate kappa_D so it works for elliptic as well as tapered (with twist). 
    kappa_D = 0.0
    for i in range(1, len(a_vector)):
        kappa_D = kappa_D + (i+1)*((a_vector[i]/a_vector[0])**2)
    print("kappa_D:\n", kappa_D, "\n")
    # print(kappa_D)

    #### now, re-evaluate kappa_L so it works for elliptic as well as tapered (with twist)
    kappa_L = (1-(1+(np.pi*wing.aspect_ratio/wing.airfoil_lift_slope))*a_vector[0])/((1+(np.pi*wing.aspect_ratio/wing.airfoil_lift_slope))*a_vector[0])

    #Now, solve for C_L_alpha of the entire wing (include the choice of elliptic vs rectangular/tapered)
    C_L_alpha = (wing.airfoil_lift_slope)/((1+(wing.airfoil_lift_slope/(np.pi*wing.aspect_ratio)))*(1 + kappa_L))
    print("C_L_alpha: \n", C_L_alpha, "\n")  


    #### now allow the calculation of the optimum max washout.
    if wing.Omega == "optimum" and wing.wing_type == "elliptic":
        Omega = 0.0
        print("Elliptic Optimum Omega\n", np.degrees(Omega))
    elif wing.Omega == "optimum" and wing.wing_type == "tapered":
        Omega = (kappa_DL*wing.CL_design)/(2*kappa_DOmega*C_L_alpha)
        print("\nOptimum Omega:", np.degrees(Omega),"\n")
    else:
        Omega = np.radians(wing.Omega)
        print("\nOmega:", np.degrees(Omega),"\n")

    
    ### now calculate the aileron distribution using straight line hinges.
    aileron_dist_array = np.zeros((wing.N_nodes))
    aileron_zb_positions = [wing.end_first_aileron_zb, wing.begin_first_aileron_zb, wing.begin_second_aileron_zb, wing.end_second_aileron_zb]
    # print("\naileron_zb_position:", aileron_zb_positions, "\n")
    aileron_thetas = np.zeros((len(aileron_zb_positions)))
    for i in range(len(aileron_thetas)):
        aileron_thetas[i] = np.arccos(2*aileron_zb_positions[i])
    # print("\naileron_thetas:", np.degrees(aileron_thetas[0]), np.degrees(aileron_thetas[1]),  np.degrees(aileron_thetas[2]),  np.degrees(aileron_thetas[3]))
    aileron_c_over_b = np.zeros((len(aileron_zb_positions)))
    if wing.wing_type == "elliptic":
        for i in range(len(aileron_thetas)):
            aileron_c_over_b[i] = ((4*1)/(np.pi*wing.aspect_ratio))*np.sin(aileron_thetas[i])
    elif wing.wing_type == "tapered":
        for i in range(len(aileron_thetas)):
            aileron_c_over_b[i] = (2*1/(wing.aspect_ratio*(1 + wing.taper_ratio)))*(1-((1-wing.taper_ratio)*abs(np.cos(aileron_thetas[i]))))
    # print("\naileraon_c_over_b:", aileron_c_over_b, "\n")

    # now, determine the cfb values. 
    cfb_array = np.zeros((len(aileron_zb_positions)))
    aileron_cfc_positions = [wing.end_aileron_cfc, wing.begin_aileron_cfc, wing.begin_aileron_cfc, wing.end_aileron_cfc]
    for i in range(len(aileron_thetas)):
        cfb_array[i] = aileron_c_over_b[i]*(-0.75 + aileron_cfc_positions[i])
    
    # print("\ncfb_array:", cfb_array, "\n")
    first_aileron_slope = (cfb_array[0]-cfb_array[1])/(aileron_zb_positions[0]-aileron_zb_positions[1])
    # print("first aileron slope:", first_aileron_slope)
    second_aileron_slope = (cfb_array[3]-cfb_array[2])/(aileron_zb_positions[3]-aileron_zb_positions[2])     
    # print("second aileron slope:", second_aileron_slope)
    

    chi_array = np.zeros((wing.N_nodes))
    cfc_array = np.zeros((wing.N_nodes))
    for i in range(0, wing.N_nodes):
        if z_b_array[i] < aileron_zb_positions[3]:
            # print("\nThis z_b:", z_b_array[i], "\n")
            chi_array[i] = 0.0
            cfc_array[i] = 0.0

        elif z_b_array[i] >= aileron_zb_positions[3] and z_b_array[i] <= aileron_zb_positions[2]:
            # print("\nThis z_b:", z_b_array[i], "\n")
            y_val = second_aileron_slope*(z_b_array[i] - aileron_zb_positions[3]) + cfb_array[3]
            cfc = (y_val/(chord_dist_array[i])) + 0.75
            cfc_array[i] = cfc

            theta_f = np.arccos(2*cfc - 1)
            ideal_flap_effectiveness = 1 - ((theta_f-np.sin(theta_f))/(np.pi))
            section_flap_effectiveness = wing.hinge_efficiency*wing.deflection_efficiency*ideal_flap_effectiveness
            chi_array[i] = section_flap_effectiveness
        
        elif z_b_array[i] >= aileron_zb_positions[2] and z_b_array[i] <= aileron_zb_positions[1]:
            chi_array[i] = 0.0
            cfc_array[i] = 0.0

        
        elif z_b_array[i] >= aileron_zb_positions[1] and z_b_array[i] <= aileron_zb_positions[0]:
            # print("\nThis z_b:", z_b_array[i], "\n")
            y_val = first_aileron_slope*(z_b_array[i] - aileron_zb_positions[0]) + cfb_array[0]
            # print("\nchord here:\n", chord_dist_array[i])
            # print("\ny_val_here:", y_val,"\n")
            cfc = (y_val/chord_dist_array[i]) + 0.75
            cfc_array[i] = cfc
            # print("\ncfc_here:", cfc, "\n")
            theta_f = np.arccos(2*cfc - 1)
            # print("\ntheta_f:\n", theta_f)
            ideal_flap_effectiveness = 1 - ((theta_f-np.sin(theta_f))/(np.pi))
            # print("\nidea flap effectiveness:\n", ideal_flap_effectiveness)
            section_flap_effectiveness = wing.hinge_efficiency*wing.deflection_efficiency*ideal_flap_effectiveness
            # print("\nsection flap effectiveness:\n", section_flap_effectiveness)
            chi_array[i] = -section_flap_effectiveness
        
        elif z_b_array[i] >aileron_zb_positions[0]:
            chi_array[i] = 0.0
    # mult = -1
    # chi_array = mult*chi_array     
    print("\nchi_array:\n", chi_array)
    print("\ncfc_array:\n", cfc_array)


    #### now, calculate the cn array by multiplying the C_inverse matrix by the values from the chi distribution
    c_vector = np.matmul(C_matrix_inv, -chi_array)
    # add the c vector to the solutions text file along with the other stuff from before

    #### now, calculate the dn array by multiplying the C_inverse matrix by the values from the following cos_theta distribution
    cos_theta_array = np.zeros((len(theta_array)))
    for i in range(len(theta_array)):
        cos_theta_array[i] = np.cos(theta_array[i])
    d_vector = np.matmul(C_matrix_inv, cos_theta_array)
    # add the d vector to the solutions text file along with the other stuff from before


    #### to calculate the rolling and yawing moments, we first need cl_da, and cl_p_bar
    Cl_da = -((np.pi*wing.aspect_ratio)/(4))*c_vector[1] # eq 1.8.60 mechanics of flight
    print("\nCl_da:", Cl_da, "\n")
    Cl_p_bar = -((np.pi*wing.aspect_ratio)/(4))*d_vector[1] # eq 1.8.61 mechanics of flight
    print("\nCl_p_bar:", Cl_p_bar, "\n")

    if wing.pbar == "steady":
        p_bar = (-Cl_da/Cl_p_bar)*wing.aileron_deflection
        print("\npbar steady:", p_bar, "\n")
    else:
        p_bar = wing.pbar
        print("\npbar:", p_bar, "\n")


    # this calculates the rolling moment coefficient
    Cl = Cl_da*wing.aileron_deflection + Cl_p_bar*p_bar
    print("\nRolling moment coefficient, Cl", Cl, "\n")


    # now find A_vector which takes into account the fourier solutions from a_n, b_n, c_n, and d_n
    A_vector = np.zeros((len(a_vector)))
    for i in range(len(A_vector)):
        A_vector[i] = a_vector[i]*(wing.alpha_root-0.0) - b_vector[i]*Omega + c_vector[i]*(wing.aileron_deflection) + d_vector[i]*(p_bar) # 1.8.50 mech of flight
    print("\nA_vector:", A_vector, "\n")


    # now calculate the yawing moment coefficient
    Cnstart = (np.pi*wing.aspect_ratio)/4
    Csubtract = ((np.pi*wing.aspect_ratio*p_bar)/8)*(A_vector[0]+A_vector[2])
    Cnsum = 0
    for i in range(1, wing.N_nodes):
        Cnsum = Cnsum + (2*(i+1)-1)*(A_vector[i-1]*A_vector[i])
    Cn = (Cnstart*Cnsum) - Csubtract # eq 6.8 aeronautics engineering handbook
    print("\nCn:", Cn, "\n")

    
    #Now, find the lift coefficient at the operating angle of attack
    C_L = C_L_alpha*((wing.alpha_root-0.0)-epsilon_omega*Omega)   ##### assuming a zero lift angle of attack due to a symmetric airfoil
    print("\nC_L: \n", C_L, "\n")
    
    #Now, find the induced drag coefficient that neglects roll and aileron
    C_Di_neglect_ail_and_roll = ((C_L**2)*(1+kappa_D)-kappa_DL*C_L*C_L_alpha*Omega+kappa_DOmega*(C_L_alpha*Omega)**2)/(np.pi*wing.aspect_ratio) # eq 1.8.25 mech of flight
    print("C_Di neglecting aileron deflection and roll: \n", C_Di_neglect_ail_and_roll, "\n") 

    #Now, find the induced drag coefficient that accounts for aileron deflection and roll. 
    CDistart = np.pi*wing.aspect_ratio
    CDisubtract = ((np.pi*wing.aspect_ratio*p_bar)/2)*(A_vector[1])
    CDsum = 0.0
    for i in range(0, wing.N_nodes):
        CDsum = CDsum + ((i+1)*A_vector[i]**2) 
    CDi_include_ail_and_roll = (CDistart*CDsum) - CDisubtract # eq 6.6 aeronautics engineering handbook
    print("C_Di including aileron deflection and roll: \n", CDi_include_ail_and_roll, "\n") 


    #export answers to text file
    with open("Solutions.txt", "w") as file:
        file.write("C_matrix:\n")
        np.savetxt(file, C_matrix, fmt='%.11f', delimiter='\t')

        file.write("\nC_matrix_inverse:\n")
        np.savetxt(file, C_matrix_inv, fmt='%.11f', delimiter='\t')

        file.write("\nFourier Coefficients (a_n):\n")
        np.savetxt(file, a_vector, fmt='%.11f', delimiter='\t')

        file.write("\nFourier Coefficients (b_n):\n")
        np.savetxt(file, b_vector, fmt='%.11f', delimiter='\t')

        file.write("\nFourier Coefficients (c_n):\n")
        np.savetxt(file, c_vector, fmt='%.11f', delimiter='\t')

        file.write("\nFourier Coefficients (d_n):\n")
        np.savetxt(file, d_vector, fmt='%.11f', delimiter='\t')

        file.write("\nFourier Coefficient Sum (A_n):\n")
        np.savetxt(file, A_vector, fmt='%.11f', delimiter='\t')

    ####now, get to plotting, 
    # first, determine the location of the ends based on the quarter-chord being centered at zero.
    leading_edge_array = np.zeros((len(chord_dist_array)))
    for i in range(0, wing.N_nodes):
        leading_edge_array[i] = 0.25*(chord_dist_array[i])
    # print("leading_edge_array:, \n", leading_edge_array, "\n")

    #now, determine the trailing edge location based on quarter chord and chord distribution
    trailing_edge_array = np.zeros((len(chord_dist_array)))
    for i in range(0, wing.N_nodes):
        trailing_edge_array[i] = -0.75*(chord_dist_array[i])
    # print("trailing_edge_array:, \n", trailing_edge_array, "\n")

    #now, explicity define the chord line x coordinates and y coordinates
    chord_x = [-0.5, 0.5]
    chord_y = [0.0, 0.0]
    
    # plot the quarter chord (lifting line)
    plt.plot(chord_x, chord_y, label = "Lifting Line")

    if wing.wing_type == "tapered":
        # Indices for the first, middle, and last values
        indices_to_plot = [0, len(z_b_array) // 2, len(z_b_array) - 1]

        # Select the values at the specified indices
        z_b_vals = [z_b_array[i] for i in indices_to_plot]
        y_leading = [leading_edge_array[i] for i in indices_to_plot]
        y_trailing = [trailing_edge_array[i] for i in indices_to_plot]

        # Create the plot for leading edge
        plt.plot(z_b_vals, y_leading, label='Leading Edge')

        # Create the plot for trailing edge
        plt.plot(z_b_vals, y_trailing, label='Trailing Edge')

        # connect the ends of the wings. 
        plt.plot([z_b_vals[0], z_b_vals[0]], [y_leading[0], y_trailing[0]], color='gray')
        plt.plot([z_b_vals[2], z_b_vals[2]], [y_leading[2], y_trailing[2]], color='gray')

    
    elif wing.wing_type == "elliptic":
        plt.plot(z_b_array, leading_edge_array, label = 'Leading Edge')
        plt.plot(z_b_array, trailing_edge_array, label = 'Trailing Edge')

    # plot the nodes along the lifting line.
    plt.plot([z_b_array[0], z_b_array[0]], [leading_edge_array[0], trailing_edge_array[0]], label = "Nodes", color = 'red')
    for i in range(1, len(z_b_array)):
        plt.plot([z_b_array[i], z_b_array[i]], [leading_edge_array[i], trailing_edge_array[i]], color = 'red')


    # plot the spanwise part of the ailerons on the planform itself. 
    plt.plot([aileron_zb_positions[0], aileron_zb_positions[1]], [cfb_array[0], cfb_array[1]], label = "Ailerons", color = "black")
    plt.plot([aileron_zb_positions[2], aileron_zb_positions[3]], [cfb_array[2], cfb_array[3]], color = "black")

    # now, plot the lines down to the trailing edge chord distribution to complete the ailerons
    for i in range(len(aileron_zb_positions)):
        plt.plot([aileron_zb_positions[i], aileron_zb_positions[i]], [cfb_array[i], -0.75*aileron_c_over_b[i]], color = "black")

    # Add labels and a legend (optional)
    plt.xlabel('z/b')
    plt.ylabel('C/b')
    plt.title('Planform')
    plt.gca().set_aspect('equal') # making the plotting axis scales equal to eachother. 
    plt.legend()
    plt.ylim(-0.5,0.5)
    plt.figure()
    # plt.show()

    #plot the washout distribution
    plt.plot(z_b_array, washout_distribution)
    plt.title('Washout Distribution')
    plt.xlabel('z/b')
    plt.ylabel('Omega')

    plt.gca().set_aspect('equal') # making the plotting axis scales equal to eachother. 
    plt.figure()

    # plot the aileron distribution chi(z/b)
    plt.plot(z_b_array, chi_array, label = "Ailerons")
    plt.title('Chi distribution')
    plt.xlabel('z/b')
    plt.ylabel('Chi')
    # plt.figure()
    plt.show()

###### Go and run the test cases plotting induced drag as a function of CL 
