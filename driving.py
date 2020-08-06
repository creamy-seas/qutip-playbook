#!/usr/bin/env python
# coding: utf-8

# In[4]:


from qutip import *
from scipy import *
import os;
from colour import Color
import numpy as np
import matplotlib.pyplot as plt


# # We will be working with frequencies, and NOT energies 
# 
# ## We can switch to density matrix with the ket2dm(ket) functions

# In[5]:


#standard states to use initially
stt_excited = basis(2,0);
stt_relaxed = basis(2,1);
stt_supp = (basis(2,0) + basis(2,1)).unit();
stt_suppm = (basis(2,0) - basis(2,1)).unit();
stt_suppi = (basis(2,0) + 1.j * basis(2,1)).unit();
stt_suppim = (basis(2,0) - 1.j * basis(2,1)).unit();
stt_rnd = (basis(2,0) + 2 * basis(2,1)).unit();


# In[6]:


def prepare_run(run_name):
    """
    Prepares Bloch sphere and subdirectories when making animations
    """
    #create sphere
    b = Bloch();
    b.vector_color = ['r'];
    b.point_color = grad_red_blue;
    b.point_size = [50];
    
    #create folder
    bashcommand = "rm -rf " + run_name;
    os.system(bashcommand);
    bashcommand = "mkdir " + run_name;
    os.system(bashcommand);
    
    return b;

def save_bloch(plot_object, directory):
    """
    Save a plot and clear it
    """
    plot_object.save(dirc = directory);
    plot_object.clear();
    
def evl_map(H, stt, blch, no_points, run_name):
    """
    Evaluate the evolution under H of the initial state stt
    
    H    Hamiltonian of the system
    stt  The initial state of the system
    blch Bloch sphere
    no_points is number of frames to make
    run_name is unique identifier
    """
    evl_map_frq(H, stt, blch, no_points, run_name, 1);

def evl_map_freq(H, stt, blch, no_points, run_name, freq):
    """
    Evaluate the evolution under H of the initial state stt

    H    Hamiltonian of the system
    stt  The initial state of the system
    blch Bloch sphere
    no_points is number of frames to make
    run_name is unique identifier
    frq  gives the frequency of rotation - can be used for a smooth animation
    """
    temp_point_x = []; temp_point_y = []; temp_point_z = [];
    for i in range(no_points):
        period = 2*pi/freq;
        time = 2*period*float(i)/no_points;
        
        #apply the evolution operator onto the state
        U = (-1.j * H * time).expm();
        evolved_state = U * stt;
        
        x = expx(evolved_state); y = expy(evolved_state); z = expz(evolved_state);
        temp_point_x.append(x);
        temp_point_y.append(y);
        temp_point_z.append(z);

        blch.add_vectors([x, y, z]);
        blch.add_points([temp_point_x, temp_point_y, temp_point_z], 'm');
        save_bloch(blch, run_name);
        
def evl_map_freq_unitary(H, stt, blch, no_points, run_name, freq, unitary_func):
    """
    Evaluate the evolution under H of the initial state stt.
    If a unitary transformation has been applied to get into a rotating frame, then 
    the state vectors are inverted!

    H    Hamiltonian of the system
    stt  The initial state of the system
    blch Bloch sphere
    no_points is number of frames to make
    run_name is unique identifier
    frq  gives the frequency of rotation - can be used for a smooth animation
    unitary_func  is the unitary transformation that we applied to get H
    """
    temp_point_x = []; temp_point_y = []; temp_point_z = [];
    for i in range(no_points):
        period = 2*pi/freq;
        time = 2*period*float(i)/no_points;
        
        #apply the evolution operator onto the state
        U = (-1.j * H * time).expm();
        evolved_state = U * stt;

        #correct for the untiary evolution we applied
        U0_dag = unitary_func(time).dag();
        evolved_state = U0_dag * evolved_state;
        
        x = real(expx(evolved_state)); y = real(expy(evolved_state)); z = real(expz(evolved_state));
        temp_point_x.append(x);
        temp_point_y.append(y);
        temp_point_z.append(z);

        blch.add_vectors([x, y, z]);
        blch.add_points([temp_point_x, temp_point_y, temp_point_z], 'm');
        save_bloch(blch, run_name);
        
def evl_master(H, stt, blch, no_points, run_name, dissapation):
    """
    Evaluate the evolution under H of the initial state stt

    H    Hamiltonian of the system
    stt  The initial state of the system
    blch Bloch sphere
    no_points is number of frames to make
    run_name is unique identifier
    dissapation is the Linbland term
    """
    times = linspace(0, 15, no_points);
    #sometimes there are more than one type of dissapation
    if(isinstance(dissapation, (list))):
        result = mesolve(H, stt, times, dissapation, [sigmax(), sigmay(), sigmaz()]);
    else:
        result = mesolve(H, stt, times, [dissapation], [sigmax(), sigmay(), sigmaz()]);
    temp_point_x = []; temp_point_y = []; temp_point_z = [];
    
    for i in range(no_points):
        
        temp_point_x.append(result.expect[0][i]);
        temp_point_y.append(result.expect[1][i]);
        temp_point_z.append(result.expect[2][i]);

        blch.add_vectors([result.expect[0][i], result.expect[1][i], result.expect[2][i]]);
        blch.add_points([temp_point_x, temp_point_y, temp_point_z], 'm');
        
        save_bloch(blch, run_name);


# In[14]:


#create colour gradient from red to blue
red = Color("red");
temp_red_blue = list(red.range_to(Color("blue"),50));
grad_red_blue = [];
for i in range(len(temp_red_blue)):
    grad_red_blue.append(temp_red_blue[i].hex_l); #extracting the hex values

#expectation values for the Bloch sphere
def expx(sst):
    return expect(sigmax(),sst);
def expy(sst):
    return expect(sigmay(),sst);
def expz(sst):
    return expect(sigmaz(),sst);

#plotting and saving
def prepare_run_detailed(run_name):
    """
    Create a bloch sphere and 2 axes to plot the sigmax and sigmaz evolution on
    """
    #create axes and figure
    fig = plt.figure(figsize=(6,10));
    ax_blch = fig.add_axes([0., 0., 2, 1.4], projection='3d'); 
    blch = Bloch(axes=ax_blch);
    ax0 = fig.add_axes([2.2, 0.75, 2, 0.5]);
    ax1 = fig.add_axes([2.2, 0.15, 2, 0.5]);

    #format the bloch sphere
    blch.vector_color = ['r'];
    blch.point_color = grad_red_blue;
    blch.point_size = [50];
    
    #format the axes
    ax1.set_xlabel('time (arb)', size = 20)
    plt.rc('xtick', labelsize = 20); 
    plt.rc('ytick', labelsize = 20); 
    plt.rc('ytick.major', size = 8, pad = 14)
    ax0.set_xticklabels([]);
    ax0.get_shared_x_axes().join(ax0, ax1);
    ax0.set_ylim([-1,1])
    ax1.set_ylim([-1,1])
    ax1.set_xlim([0,1])
    ax0.set_xlim([0,1])
    for i in [ax0, ax1]:
        i.set_yticks([-1, -0.5, 0, 0.5, 1])
        i.set_yticklabels(["-1", "-0.5", "0", "0.5", "1"]) 

    #global param
    ax0.set_title("$\sigma_z$", size = 40, pad = 20, color = 'red');
    ax1.set_title("$\sigma_x$", size = 40, pad = 20, color = 'blue');
    ax0.grid(True, which='both');
    ax1.grid(True);

    #create folder
    bashcommand = "rm -rf " + run_name;
    os.system(bashcommand);
    bashcommand = "mkdir " + run_name;
    os.system(bashcommand);
    
    return [blch, fig];

def save_bloch_detailed(folder, index, figure, blch):
    location = folder + "/bloch_" + str(index) + ".png";
    figure.savefig(location, bbox_inches = 'tight');
    
def evl_detailed(H, stt, blch, no_points, run_name, freq, fig, unitary):
    """
    Evaluate the evolution under H of the initial state stt 

    H    Hamiltonian of the system
    stt  The initial state of the system
    blch Bloch sphere
    no_points is number of frames to make
    run_name is unique identifier
    frq  gives the frequency of rotation - can be used for a smooth animation
    fig  is the figure with the x (1) and z (2) axes
    """
    
    ax = fig.axes;
    temp_point_x = []; temp_point_y = []; temp_point_z = [];
    
    for i in range(no_points):
        period = 2*pi/freq;
        time = 2*period*float(i)/no_points;
        
        #apply the evolution operator onto the state
        U = (-1.j * H * time).expm();
        evolved_state = U * stt;
        
        #find the expectation values
        x = real(expx(evolved_state)); y = real(expy(evolved_state)); z = real(expz(evolved_state));
        temp_point_x.append(x);
        temp_point_y.append(y);
        temp_point_z.append(z);

        #plot on the bloch sphere
        blch.clear();
        blch.add_vectors([float(x), float(y), float(z)]);
        blch.add_points([temp_point_x, temp_point_y, temp_point_z], 'm');
        blch.make_sphere();
        
        #plot on the axes
        ax[1].plot(float(1)/no_points*i, z, 'o', c = 'red');
        ax[2].plot(float(1)/no_points*i, x, 'o', c = 'blue');
        
        save_bloch_detailed(run_name, i, fig, blch);
        
def evl_detailed_unitary(H, stt, blch, no_points, run_name, freq, unitary_func):
    """
    Evaluate the evolution under H of the initial state stt.
    If a unitary transformation has been applied to get into a rotating frame, then 
    the state vectors are inverted!

    H    Hamiltonian of the system
    stt  The initial state of the system
    blch Bloch sphere
    no_points is number of frames to make
    run_name is unique identifier
    frq  gives the frequency of rotation - can be used for a smooth animation
    fig  is the figure with the x (1) and z (2) axes
    unitary_func  is the unitary transformation that we applied to get H
    """
    ax = fig.axes;
    temp_point_x = []; temp_point_y = []; temp_point_z = [];
    
    for i in range(no_points):
        period = 2*pi/freq;
        time = 2*period*float(i)/no_points;
        
        #apply the evolution operator onto the state
        U = (-1.j * H * time).expm();
        evolved_state = U * stt;
        #correct for the untiary evolution we applied
        U0_dag = unitary_func(time).dag();
        evolved_state = U0_dag * evolved_state;
        
        #find the expectation values
        x = real(expx(evolved_state)); y = real(expy(evolved_state)); z = real(expz(evolved_state));
        temp_point_x.append(x);
        temp_point_y.append(y);
        temp_point_z.append(z);

        #plot on the bloch sphere
        blch.clear();
        blch.add_vectors([float(x), float(y), float(z)]);
        blch.add_points([temp_point_x, temp_point_y, temp_point_z], 'm');
        blch.make_sphere();
        
        #plot on the axes
        ax[1].plot(float(1)/no_points*i, z, 'o', c = 'red');
        ax[2].plot(float(1)/no_points*i, x, 'o', c = 'blue');
        
        save_bloch_detailed(run_name, i, fig, blch);
        
        
def evl_master_detailed(H, stt, blch, no_points, run_name, dissapation, fig):
    """
    Evaluate the evolution under H of the initial state stt and plot the sigmaz and sigmax alongside

    H    Hamiltonian of the system
    stt  The initial state of the system
    blch Bloch sphere
    no_points is number of frames to make
    run_name is unique identifier
    dissapation is the Linbland term
    fig  is the figure with the x (1) and z (2) axes
    """
    times = linspace(0, 15, no_points);
    
    #sometimes there are more than one type of dissapation
    if(isinstance(dissapation, (list))):
        result = mesolve(H, stt, times, dissapation, [sigmax(), sigmay(), sigmaz()]);
    else:
        result = mesolve(H, stt, times, [dissapation], [sigmax(), sigmay(), sigmaz()]);
    
    ax = fig.axes;
    temp_point_x = []; temp_point_y = []; temp_point_z = []; time_point = [];
    
    
    for i in range(no_points):    
        temp_point_x.append(result.expect[0][i]);
        temp_point_y.append(result.expect[1][i]);
        temp_point_z.append(result.expect[2][i]);
        
        #plot on the bloch sphere
        blch.clear();
        blch.add_vectors([result.expect[0][i], result.expect[1][i], result.expect[2][i]]);
        blch.add_points([temp_point_x, temp_point_y, temp_point_z], 'm');
        blch.make_sphere();
        
        #plot on the axes
        ax[1].plot(float(1)/no_points*i, result.expect[2][i], 'o', c = 'red');
        ax[2].plot(float(1)/no_points*i, result.expect[0][i], 'o', c = 'blue');
        
        save_bloch_detailed(run_name, i, fig, blch);
        
def evl_animate(run_name):
    """
    Makes and animation from the images in the run_name subdirectory
    """
    bashcommand = "rm " + run_name +".mp4";
    os.system(bashcommand);
    bashcommand = "ffmpeg -r 10 -i " + run_name + "/bloch_%01d.png " + run_name +".mp4";
    os.system(bashcommand); 
    bashcommand = "rm -r " + run_name;


# ## No dissapation, no drive

# In[77]:


### Dynamics with no disspation no drive (we will be uisng linear vectors)
no_points = 40;
run_list = [];

run_list.append(["noDrive_noDissapataion", ket2dm(stt_rnd)]);

#system state and Hamiltonian
omega = 2;
H = float(omega)/2*sigmaz();

#create animation
counter = 0;
for run in run_list:
    #extract run elements
    run_name = run[0];
    initial_state = run[1];
    
    #prepare sphere and map evolution
    blch, fig = prepare_run_detailed(run_name);
    evl_detailed(H, initial_state, blch, no_points, run_name, float(omega), fig);
    evl_animate(run_name);
    
    counter = counter + 1;
    print("Finished run " + str(counter) + "/" + str(len(run_list)));


# # No dissapation, drive

# In[79]:


## Dynamics with no disspation but with drive
no_points = 40;
run_list = [];

run_list.append(["drive_noDissapation", stt_relaxed]);

#system state and Hamiltonian. 
#here we enter the rotating frame at the frequency of the drive (energy splitting of the atom) 
#and remove the fast rotating terms
#now we are considering the system where |0> and |1> are at the same energy, and their superpositions
#are the new eigenstates of the system
omega = float(2);
Omega = float(3);
H = Omega/2*sigmax();

#create animation
counter = 0;
for run in run_list:
    #extract run elements
    run_name = run[0];
    initial_state = run[1];
    
    #prepare sphere and map evolution
    blch, fig = prepare_run_detailed(run_name);
    evl_detailed(H, initial_state, blch, no_points, run_name, float(omega), fig);
    evl_animate(run_name);
    
    counter = counter + 1;
    print("Finished run " + str(counter) + "/" + str(len(run_list)));


# In[ ]:


### Same as above, but we exit the rotating frame
no_points = 40;

run_list = [];
run_list.append(["drive_noDissapation_ROTATING", stt_relaxed]);

#system state and Hamiltonian. 
#here we enter the rotating frame at the frequency of the drive (== energy splitting of the atom) 
#and remove the fast rotating terms
#now we are considering the system where |0> and |1> are at the same energy, and their superpositions
#are the new eigenstates of the system
initial_state = stt_relaxed;
omega = float(2);
Omega = float(10);
H = Omega/2*sigmax();
def unitary_00(t):
    #This unitary transformation put us into the rotating frame , in which |0> and |1> are at the same energies
    return (-1.0j * 2 / 2 * t * sigmaz()).expm();

#prepare sphere and map evolution
blch = evl_prepare(run_name); 
evl_map_freq_unitary(H, initial_state, blch, no_points, run_name, 4, unitary_00)
evl_animate(run_name)


# ## Dissapation, No drive

# In[28]:


no_points = 60;
run_list = [];

#run_list.append(["noDrive_@Superposition_Relaxation_Dephasing", ket2dm(stt_supp), 0.3, 0.1]);
#run_list.append(["noDrive_@Superposition_Relaxation", ket2dm(stt_supp), 0.1, 0]);
#run_list.append(["noDrive_@Superposition_Dephasing", ket2dm(stt_supp), 0, 0.1]);
#run_list.append(["noDrive_@Excited_Relaxation_Dephasing", ket2dm(stt_excited), 0.1, 0.1]);
#run_list.append(["noDrive_@Excited_Relaxation", ket2dm(stt_excited), 0.1, 0]);

#raw Hamiltonian of the two level system
omega = 2;
H = float(omega)/2 * sigmaz();

#repeat for all runs
counter = 0;
for run in run_list:
    #extract run elements
    run_name = run[0];
    initial_state = run[1];
    Gamma_1 = run[2];
    Gamma_Phi = run[3];
    
    #define dissapation
    dissapation = [np.sqrt(Gamma_1)*sigmam(), np.sqrt(Gamma_Phi/2)*sigmaz()];

    #prepare sphere and map evolution
    blch, fig = prepare_run_detailed(run_name);
    evl_master_detailed(H, initial_state, blch, no_points, run_name, dissapation, fig);
    evl_animate(run_name);
    
    counter = counter + 1;
    print("Finished run " + str(counter) + "/" + str(len(run_list)));


# ## Dissapation, No drive, ROTATING FRAME

# In[16]:


no_points = 10;
run_list = [];
run_list.append(["test", ket2dm(stt_supp), 0.3, 0.1]);
#run_list.append(["noDrive_@Superposition_Relaxation_Dephasing_ROTATING", ket2dm(stt_supp), 0.3, 0.1]);
#run_list.append(["noDrive_@Superposition_Relaxation_ROTATING", ket2dm(stt_supp), 0.3, 0]);
#run_list.append(["noDrive_@Superposition_Dephasing_ROTATING", ket2dm(stt_supp), 0, 0.1]);
#run_list.append(["noDrive_@Excited_Relaxation_Dephasing_ROTATING", ket2dm(stt_excited), 0.1, 0.1]);
#run_list.append(["noDrive_@Excited_Relaxation_ROTATING", ket2dm(stt_excited), 0.1, 0]);

#raw Hamiltonian of the two level system
omega = 2;
H = float(omega)/2 * sigmaz() - omega/2 * sigmaz(); #in the rotating frame, we negatve the effect of this Hamiltonian
#repeat for all runs
counter = 0;

for run in run_list:
    #extract run elements
    run_name = run[0];
    initial_state = run[1];
    Gamma_1 = run[2];
    print(run[2])
    Gamma_Phi = run[3];
    
    #define dissapation
    dissapation = [np.sqrt(Gamma_1)*sigmam(), np.sqrt(Gamma_Phi/2)*sigmaz()];

    #prepare sphere and map evolution
    blch, fig = prepare_run_detailed(run_name);
    evl_master_detailed(H, initial_state, blch, no_points, run_name, dissapation, fig);
    evl_animate(run_name);
    
    counter = counter + 1;
    print("Finished run " + str(counter) + "/" + str(len(run_list)));


# ## Dissapation, Driving

# In[42]:


no_points = 3;
run_list = [];

run_list.append(["drive_@Excited_Relaxation_Dephasing", ket2dm(stt_excited), 0.2, 0.1]);
"""run_list.append(["drive_@Excited_Relaxation", ket2dm(stt_excited), 0.2, 0]);
run_list.append(["drive_@Excited_Dephasing", ket2dm(stt_excited), 0, 0.4]);
run_list.append(["drive_@Superposition_Relaxation_Dephasing", ket2dm(stt_supp), 0.6, 0.1]);
run_list.append(["drive_@Superposition_Relaxation", ket2dm(stt_supp), 0.3, 0]);
run_list.append(["drive_@Superposition_Dephasing", ket2dm(stt_supp), 0, 0.2]);
"""

#here we enter the rotating frame at the frequency of the drive (energy splitting of the atom) 
#and remove the fast rotating terms
omega = 2;
Omega = 1.5;
phi = 0;
H = float(Omega)/2*(sin(phi)*sigmay() + cos(phi)*sigmax());

#repeat for all runs
counter = 0;
for run in run_list:
    #extract run elements
    run_name = run[0];
    initial_state = run[1];
    Gamma_1 = run[2];
    Gamma_Phi = run[3];
    
    #define dissapation (check Ilya Notes summary if not sure)
    dissapation = [np.sqrt(Gamma_1)*sigmam(), np.sqrt(Gamma_Phi/2)*sigmaz()];

    #prepare sphere and map evolution
    blch, fig = prepare_run_detailed(run_name);
    evl_master_detailed(H, initial_state, blch, no_points, run_name, dissapation, fig);
    evl_animate(run_name);
    
    counter = counter + 1;
    print("Finished run " + str(counter) + "/" + str(len(run_list)));


# In[2]:





# In[ ]:




