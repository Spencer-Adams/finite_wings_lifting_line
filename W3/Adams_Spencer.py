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
            self.N_nodes = (2* self.nodes_per_semispan) - 1
            # self.alphas = input_vals["condition"]["alpha_root[deg]"]
            self.alpha_root = np.radians(input_vals["condition"]["alpha_root[deg]"])
            self.planform = input_vals["view"]["planform"]
            self.is_washout_distribution = input_vals["view"]["washout_distribution"]


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
        theta_array[i] = ((i)*np.pi)/(wing.N_nodes-1)
    # print("\n") 
    # print("theta")  
    # print(theta_array)
    # print("\n")

    # calculate array of z/b    
    z_b_array = np.zeros((wing.N_nodes))
    z_b_array[0] = 0.5
    z_b_array[wing.N_nodes - 1] = -0.5
    # now, calculate the middle z values
    for i in range(1, wing.N_nodes -1):
        z_b_array[i] = (1/2)*np.cos(theta_array[i]) # 1/2 because b = 1
    # print("z_b_array")
    # print(z_b_array)
    # print("\n")

    # calculate chord distribution array
    chord_dist_array = np.zeros((wing.N_nodes))
    for i in range(0, wing.N_nodes):
        if wing.wing_type == "elliptic":
            chord_dist_array[i] = ((4*1)/(np.pi*wing.aspect_ratio))*np.sin(theta_array[i])
        elif wing.wing_type == "tapered":
            chord_dist_array[i] = (2*1/(wing.aspect_ratio*(1 + wing.taper_ratio)))*(1-((1-wing.taper_ratio)*abs(np.cos(theta_array[i]))))
    # print("chord_distribution")
    # print(chord_dist_array)
    # print("\n")

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


    #### now allow the calculation of the optimum washout.
    if wing.Omega == "optimum" and wing.wing_type == "elliptic":
        Omega = 0.0
        print("Elliptic Optimum Omega\n", np.degrees(Omega))
    elif wing.Omega == "optimum" and wing.wing_type == "tapered":
        Omega = (kappa_DL*wing.CL_design)/(2*kappa_DOmega*C_L_alpha)
        print("\nOptimum Omega:", np.degrees(Omega),"\n")
    else:
        Omega = np.radians(wing.Omega)
        print("\nOmega:", np.degrees(Omega),"\n")

    #Now, find the lift coefficient at the operating angle of attack
    C_L = C_L_alpha*((wing.alpha_root-0.0)-epsilon_omega*Omega)   ##### assuming a zero lift angle of attack due to a symmetric airfoil
    print("\nC_L: \n", C_L, "\n")
    
    #Now, find the induced drag coefficient
    C_Di = ((C_L**2)*(1+kappa_D)-kappa_DL*C_L*C_L_alpha*Omega+kappa_DOmega*(C_L_alpha*Omega)**2)/(np.pi*wing.aspect_ratio)
    print("C_Di: \n", C_Di, "\n") 

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

    

    #now, get to plotting, 
    
        


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

    # Add labels and a legend (optional)
    plt.xlabel('z/b')
    plt.ylabel('C/b')
    plt.title('Planform')
    plt.gca().set_aspect('equal') # making the plotting axis scales equal to eachother. 
    plt.legend()
    plt.ylim(-0.5,0.5)
    plt.figure()
    # plt.show()

    plt.plot(z_b_array, washout_distribution)
    plt.gca().set_aspect('equal') # making the plotting axis scales equal to eachother. 
    plt.show()

###### Go and run the test cases plotting induced drag as a function of CL 
