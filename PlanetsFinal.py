"""
A code to simulate the orbits of the planets in the solar system, the Sun, the Moon and Halley's comet.
It draws on two input files:
1 -  holds timestep, length and gravity constant.
2 - holds the planets' mass, initial positions and initial velocities.

Authors: Lyndsey Scott and Flora Blake Parsons.
January - March 2018
"""

#Import required libraries and codes.
import numpy as np
import matplotlib.pyplot as plt
import copy #This is used to hold a static copy of position for calculating orbital periods and apo/periapses.

"""
The simulation class contains all the functions that are connected to and are required in our simulation of the solar system.
It initialises the time, the position list, planet list, the energies, the centre of mass of the system and an output file for the VMD simulation.
Then, it reads in the file that contains all the data for the planet list - the planets' mass, position (x,y,z) and velocity (vx,vy,vz).
It sets each of the elements of the list as the respective value, and multiplies by 1000 to fix the units. 
A function is created to correct for the centre of mass movement to stop the whole system drifting outwards due to the wobble of the Barycentre.
The positions at every step is stored in a format readable by VMD and then these are written into two files, the same but .xyz for VMD and .txt to check.
The total energy, kinetic energy and potential energy is plotted.
For each body, the orbital period is calculated.
The simulation is run. 
The force is definied as the gravitational force between two planets/objects.
The potential is defined as the gravitational potential.
"""


class Simulate: 
    def __init__(self, run_time, G): 
        self.file = str()
        self.run_time = run_time  									# How long we will run the simulation for.
        self.G = G                									#Gravity constant
        self.planet_list = []  										# Set up the planet list that is filled in make_planets.
        self.positions_list = []									# Set up the list that holds the positions of the planets. 
        self.potential_e = np.zeros(run_time, dtype=float) 			#Initialises potential energy as a numpy array of zeros.
        self.kinetic_e = np.zeros(run_time, dtype=float) 			#Initialises kinetic energy as a numpy array of zeros. 
        self.system_com = 0 #										Initialises centre of mass at the centre of the system.
        self.out_file = 'planets_out.xyz' 							#Initialises an output file for the VMD software to read. 

    def readfile(self, filein='PlanetDataFinal.txt'): 				#Reads in the file with the planets' data (their x,y,z coordinates, mass, vx,vy,vz) 
        self.file = np.loadtxt(filein)
        return self.file

    def make_planets(self):		
        """
        Stores each line of the input file as a Planet object.
        :return: List of Planet objects, as read from in file.
        """						
        self.planet_list = np.zeros(len(self.file), dtype=object) #Initialises the planet list as a numpy array of zeros, the length of the planet data file.

        for i in range(len(self.file)): 		
            mass = self.file[i][0]
            x = self.file[i][1] * 1000		#Multiplied by 1000 to change from km to m.
            y = self.file[i][2] * 1000
            z = self.file[i][3] * 1000
            vx = self.file[i][4] * 1000		#Multiplied by 1000 to change from km/s to m/s
            vy = self.file[i][5] * 1000
            vz = self.file[i][6] * 1000
            label = i  # self.file[i][7]	#Labels the objects 0,1,2,3... etc 

            self.planet_list[i] = Planet(x, y, z, vx, vy, vz, mass, label)
        return self.planet_list
        
        
    def com_correction(self):
        """
        Calculates the center of mass velocity correction, which is needed as otherwise everything will drift outwards due to the wobble of the Barycentre.
        :return: Array velocity correction vector.
        """
        p_total = 0 						#Initialises the momentum as zero
        tot_mass = 0 						#Initialises total mass as zero 
        com = 0 							#Initialises centre of mass as the origin
        for planet in self.planet_list:
            p_total += planet.momentum() 	# p_total += planet.momentum() is equivalent to p_total = p_total + planet.momentum()
            tot_mass += planet.mass  		# Allows for variable mass planets if that's what you would want --> Halley's comet will lose mass over time.
            com += np.multiply(planet.mass, planet.position)	#Centre of mass calculated as a numpy array of mass x position.
        self.system_com = np.divide(com, tot_mass)  			# Sets COM of system
        v_com = np.divide(p_total, tot_mass) 					#This is the correction - the momentum divided by the mass to give the velocity in an array. 
        for planet in self.planet_list:
            planet.velocity -= v_com			
        return v_com
        
    def store_positions(self, i):
        """
        Stores the positions of all objects in self.planet_list at a given time step
        in a format readable by VMD.
        :param i: Time step (int).
        :return: List of positions in VMD format up to current time step.
        """
        self.positions_list.append(len(self.planet_list))  #Appends the position list for each planet in the planet list.
        self.positions_list.append('Point = %d' % (i+1)) # So as to start at 1 not 0
        # print(len(self.planet_list))
        for planet in self.planet_list:
            self.positions_list.append(str(planet))
            # print(planet)
        return self.positions_list

    def write_out_file(self, out_file='planet_out.xyz'):
        """
        Writes positions list to an out file, which will then be used as the input for VMD. Default 'planet_out.xyz'.
        :param out_file: Out put files name.
        :return: 0.
        """
        file = open(out_file, 'w')
        txtfile = open('planet_out.txt', 'w')
        for item in self.positions_list:
            file.write("%s\n" % item) #%s makes it a string; \n makes a new line 
            txtfile.write("%s\n" % item)  # Text file to check output
        file.close()
        txtfile.close()
        return 0 

    def energy_plot(self):
        ek = self.kinetic_e #the kinetic energy of the planet itself
        ep = self.potential_e #potential energy of the planet itself
        e_tot = [ek[i] + ep[i] for i in range(len(ek))] #total energy is the kinetic energy + the potential energy for each timestep.
        t = [i for i in range(len(ek))] #time is a list where each increment is a timestep.
        plt.plot(t, ek, linestyle='--', color='red') #Plot
        plt.plot(t, ep, linestyle=':', color='blue')
        plt.plot(t, e_tot, color='black')
        plt.title('Energy vs Time')
        plt.legend(['Kinetic Energy', 'Potential Energy', 'Total Energy'])
        plt.xlabel('Time')
        plt.ylabel('Energy ($J$)')
        plt.show()

    #def orbital_periods(self): #calculates the orbital period of each planet using the centre of mass correction and their input information
        #for planet in self.planet_list: #calculates for each planet in the list
            #planet.calc_orbital_period(self.system_com, self.planet_list[0].mass)  If you want to compare to Kepler's law

    def run_sim(self, steps_per_out=5, dt=87600, vmd_out=False, twoD_plot=True):
        """
        Runs orbital simulation by calculating the net force on each Planet
        object then updates each planet resulting from the applied force.
        :param steps_per_out: The number of time steps between each segment is sent to an out file.
        :param dt: Time step in seconds per unit time of simulation.
        :param vmd_out: True if a VMD appropriate out file is needed.
        :param twoD_plot: If a 2D x, y plot is wanted.
        :return: 0 on completion.
        """
        run_time = self.run_time					# How long we will run for as set in the input file
        time = 0  									# Set initial time to 0
        self.com_correction()  						# Centre of mass correction.
        while time < run_time:  					# Run while current time is less than total run time, or else would be an error. 
            for i in range(len(self.planet_list)):  # running for each planet in the list
                force = 0 							# setting the initial conditions to zero
                potential_e = 0
                for j in range(len(self.planet_list)):
                    if i != j:  					# Sums force and potential energy from all other planets except itself.
                        force += self.find_force(self.planet_list[i], self.planet_list[j]) #Applies the force from the two planets and adds to the existing force.
                        potential_e += self.find_potential(self.planet_list[i], self.planet_list[j]) #Applies the potential from the two planets and adds to the existing potential.
                self.planet_list[i].force = force 
                self.potential_e[time] = potential_e
            ek = 0  								# set initial kinetic energy to zero
            for planet in self.planet_list:
                planet.update(dt=dt, sys_com=self.system_com, time=time)  # Update the positions, velocities and accelerations of each planet
                ek += planet.kinetic_energy() 							  # Update kinetic energy
                if time % 1000 == 0 and twoD_plot:						  # % operator means if time divided by 1000 has a remainder, would not have an integer amount of timesteps, which we dont want, condition Can be removed if a 2d plot is not wanted.
                    plt.plot(planet.position[0], planet.position[1], marker='o') 
            self.kinetic_e[time] = ek
            if time % steps_per_out == 0 and vmd_out:
                self.store_positions(i=time)
            time += 1  # increase time
        if twoD_plot:
            plt.show() #Plot the 2d plot if wanted
        if vmd_out:  # On completion if a VMD out file is desired then make the file.
            self.write_out_file()
        return 0 #close function

    def find_potential(self, planet, other_planet): 
        r = np.subtract(other_planet.position, planet.position) #the vector distance between the two planets
        r_mag = np.linalg.norm(r) # the scalar magnitude of the distance between the two planets
        ep = -self.G*planet.mass*other_planet.mass/r_mag #gravitational potential energy of the two planets (scalar) = Gmm/r
        return ep

    def find_force(self, planet, other_planet): 
        """
        Calculate the gravitational force between two Planet objects.
        :param planet: Planet the force is being applied on.
        :param other_planet: Planet the force is being applied by.
        :return: Gravitational force vector applied on planet.
        """
        r = np.subtract(other_planet.position, planet.position) #the vector distance between the two planets
        r_mag = np.linalg.norm(r) #scalar distance
        r_hat = np.divide(r, r_mag) #vector component of r (ie the direction)
        force_mag = self.G*planet.mass*other_planet.mass/r_mag**2   #magnitude of the gravitational force = Gmm/r^2
        force_vector = np.multiply(force_mag, r_hat) 				#makes it a vector
        force = force_vector
        return force
        

"""
The planet class contains all the functions and that are connected to and are required by the planets themselves. 
It initialises the properties of the planets: the position, previous position, acceleration and velocities as numpy arrays; the force, the angle and the orbit number as zero; the mass and seperation increase and stabilising as True.
It creates a position string suitable for VMD.
The apo/periapses are calculated and printed.
Kinetic energy, momentum and acceleration are calculated.
The Verlet integration method is used to find the next positions and velocities of the planets.
The orbital period is found using trig.
These values are updated. 
"""


class Planet(object):
    def __init__(self, x, y, z, vx, vy, vz, mass, label):  	# this initialises the properties we will go on to use
        self.label = label
        self.position = np.array([x, y, z], float) 			#Creates position as an array of the x,y and z scalers. 
        self.velocity = np.array([vx, vy, vz], float)		#Creates velocity as an array of the x,y and z speeds. 
        self.acceleration = np.array([0, 0, 0], float)		#Initialises acceleration as an array of zeros to be filled later.
        self.force = np.array([0, 0, 0], float)				#Initialises the force as an array of zeros to be filled later.
        self.mass = float(mass)								#Sets the mass as a float
        self.prev_position = np.array([x, y, z], float)		#Creates a previous position array for calculating the periods.
        self.r_increasing = True              				#Sets the seperation increase to positive for the beginning of def apo_peri
        self.theta = 0										#Initialises the angle and the orbit number as zero
        self.orb_no = 0
        self.stabilising = True								#Sets stabilisation as true for the calculation of periods and apo/periapses
        
        
    def __str__(self):
        """
        Creates a position string, suitable for a VMD output.
        :return: Position string.
        """
        string = ("%s %.2f %.2f %.2f" % (self.label, self.position[0], self.position[1], self.position[2])) #Formats the string for VMD.
        return string

    def apo_peri(self, com, time):  # The stabilisation is to counteract the effect of the wobble in centre of mass.
        """
        Prints radius of apoaps and periaps when a planet reaches it.
        :param com: Center of mass of the overall planet system.
        :return: Number
        """
        if time != 0:  # This and stabilising = true prevents planets which have radius that is decreasing initially to say it is at a periapsis.
            r_prev = np.linalg.norm(np.subtract(self.prev_position, com))   # The seperation between the centre of mass and the previous position.
            r = np.linalg.norm(np.subtract(self.position, com))				# The seperation between the centre of mass and the current position.
            if self.r_increasing:
                if r < r_prev and self.stabilising: 		#If the seperation is less than the previous seperation but it self stabilises it ISNT a periapsis hence False
                    self.r_increasing = False
                if r < r_prev and not self.stabilising:		#If the seperation is less than the previous seperation and it doesn't stabilise it IS a periapsis.
                    self.r_increasing = False
                    print("%s, at periapsis at a radius of "		#Prints the seperation.
                          "%.0fm from the system COM." % (self.label, r))
            elif not self.r_increasing and not self.stabilising:		#r_increasing is now false
                if r > r_prev:
                    self.r_increasing = True						#If the seperation is greater than the previous seperation it is a periapses.
                    print("%s, at apoapsis at a radius of "
                          "%.0fm from the system COM." % (self.label, r))
            if self.stabilising:
                self.stabilising = False

    def kinetic_energy(self):			#KE = 1/2 mv^2
        vel_mag = np.linalg.norm(self.velocity)
        ek = 0.5*self.mass*vel_mag**2
        return ek

    def momentum(self):					# p = mv
        p = np.multiply(self.mass, self.velocity)
        return p

    def calc_acceleration(self):		# F=ma
        self.acceleration = np.divide(self.force, self.mass)
        return self.acceleration

    def next_position(self, dt):		#The velocity verlet integration method for finding next postion
    	self.prev_position = self.position.copy()			#Stores a copy of the previous position that then won't change for use in the period calculation.
    	self.position += self.velocity * dt + 0.5 * (dt ** 2) * self.force / self.mass
    	return self.position

    def next_velocity(self, dt):		#The velocity verlet integration method for finding next velocity.
        self.velocity = np.add(self.velocity, self.acceleration*dt)
        return self.velocity
	   	
    def orbital_period(self, com, time, dt):		#Calculate the angle between the position and the previous position from the centre of mass and sum until they equal 2pi.
        a = np.subtract(self.position, com)			#Find the vector from the current position and the centre of mass.
        b = np.subtract(self.prev_position, com)	#The vector from the previous position and the centre of mass.
        a_mag = np.linalg.norm(a)					#The magnitude of the above vectors.
        b_mag = np.linalg.norm(b)
        norm_a_dot_b = np.dot(a, b) / np.multiply(a_mag, b_mag)		#Using the dot product: cosTHETA = (a.b)/(mag_a*mag_b)
        theta = np.arccos(norm_a_dot_b)
        self.theta += theta							#Sum to the previous theta total
        if self.theta >= np.pi*2:					#When it equals 2pi --> 
            self.orb_no += 1							#Add one to the orbit number
            orb_time = time*dt / self.orb_no			#Calculate the time elapsed (seconds)
            orb_t_years = orb_time / (3600*24*365)		#Convert to years
            print("The orbital period of %s is %.3f years." % (self.label, orb_t_years))	#Print the period
            self.theta = 0

    def update(self, dt, sys_com, time):
        """
        Updates a planet object being acted on by a force to find
        its new acceleration, position and velocity.
        :param dt: Time step.
        :param sys_com: Position of the total systems center of mass
        :return: 0
        """
        acceleration = self.calc_acceleration()
        position = self.next_position(dt)
        velocity = self.next_velocity(dt)
        self.apo_peri(com=sys_com, time=time)
        self.orbital_period(sys_com, time, dt)
        return 0

"""
This following function found the orbital period using Keplar's laws which we wanted to avoid
but was used to check our trig version was working.
    
    def calc_orbital_period(self, sys_com, central_mass, printer=True):
        mean_radius = np.linalg.norm(np.subtract(self.position, sys_com))
        numerator = mean_radius**3*4*np.pi**2
        denominator = 6.67*10**(-11)*central_mass
        period = (numerator / denominator)**(1/2)
        period_years = period / (3600*24*365)
        if printer:
            print("The orbital period of %s is %.2f years." % (self.label, period_years))
        return period
"""



"""
Main
Loads the input, and sets run time to the first element; the timestep as the second and the gravity constant as the third.
Runs the simulation: reads the planets file, makes the planet array, finds the orbital periods, runs the simulation and plots the energy.

"""


def main():
    sys_input = np.loadtxt('input.txt')
    run_time = sys_input[0]
    dt = sys_input[1]
    G = sys_input[2]
    sim = Simulate(run_time=int(run_time), G=G)  
    sim.readfile()
    sim.make_planets()
    #sim.orbital_periods() #If you want to compare to Kepler's
    sim.run_sim(vmd_out=True, dt=dt)  	#Set vmd_out = False if VMD is not required.
    sim.energy_plot()					#Can remove if energy plot not wanted. 

main() #Runs the code

        
        
        
        
        
        
        

   
