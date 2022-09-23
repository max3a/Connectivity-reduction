# Connectivity-reduction
C++ code. The program is operated under windows with a gnu complier, and RandNum.h (from http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html) and run.h provided whitin this repository to run the program.

The model circuit is composed of 2000 excitatory neurons and 400 inhibitory neurons, and probability 0.2, the neurons were randomly connected. In inhibitory population, 200 neurons are picked to receive inhibitory pulse to mimic CA1 neurons rececving input from MSPV neuron. The neurons are labelled as:

0-1999 excitatory neurons;
2000-2199 inhibitory neurons;
2200-2399 inhibitory neurons that may receive inhibitory pulse;

The model used in this work is composed of conductance-based integrate-and-fire neurons, short-term plasticity and long-term plasticity, and the parameters used are biologically plausible, justified in previous works (Brunel and Wang, 2003; Dayan and Abbott, 2005; Hempel et al., 2000; Markram et al., 1998; Mongillo et al., 2008; Pfister and Gerstner, 2006). In particular, each layer has 2000 excitatory neurons, and 400 inhibitory neurons, n_I. Each neuron in L1 and L2 received background input from 400 independent Poisson trains. 

In main.cpp, you can decide

the time interval between two recordings, t_record.
Able/disable feedforward E->I plasticity, LTP.
the theta, and gamma frequency and number of pulses applied within one theta period during rescue when feedforward E->I plasticity is not disabled, theta, gamma, pulse_number.
the numbering of the simulation, count
To run the simulation, you can run cmd in window, go to the directory folder and type mingw32-make. It will generate an exe, named "main.exe".

By running "main.exe", 9 txt files will be generated, which are

E=3.00I=8.00P=0.20_spike_{count}_0_0.00.txt (with feedforward E->I plasticity) or E=3.00I=8.00P=0.20_spike_{count}_0_0.50.txt (without feedforward E->I plasticity)
The txt has 2 columns, which are the spike time and neuron label (Each neuron has its unique label). i.e. 0-1999 are L1 excitatory neurons; 2000-2399 are L1 inhibitory neurons; 2800-4799 are L2 excitatory neurons; 4800-4999 are L2 Re inh neurons; 5000-5199 are FF inh neurons.

text_{count}_RE_AMPA.txt
The txt has 2401 colums. The 1st column is running time of simulation in unit ms. The 2-2401 columns are the recurrent AMPA current recived by neurons labelled with (column number+1998). i.e. column number 2-2001 are L2 excitatory neurons (Neuron labels are 2800-4799); column number 2002-2201 are L2 Re inh neurons (Neuron labels are 4800-4999); column number 2202-2401 are FF inh neurons (Neuron labels are 5000-5199).

text_{count}_RE_NMDA.txt
Same as text_{count}_RE_AMPA.txt, but the 2-2401 columns are recurrent NMDA current.

text_{count}_FF_AMPA.txt
Same as text_{count}_RE_AMPA.txt, but the 2-2401 columns are feedforward NMDA current.

text_{count}_FF_NMDA.txt
Same as text_{count}_RE_AMPA.txt, but the 2-2401 columns are feedforward NMDA current.

text_{count}_FF_GABA.txt
Same as text_{count}_RE_AMPA.txt, but the 2-2401 columns are recurrent GABA current.

text_{count}_V.txt
Same as text_{count}_RE_AMPA.txt, but the 2-2401 columns are membrane potential of neurons.

E=3.00I=8.00P=0.20_W_{count}_0_0.00.txt (with feedforward E->I plasticity) or E=3.00I=8.00P=0.20_W_{count}_0_0.50.txt (without feedforward E->I plasticity)
A matrix representing the final synaptic weights after learning. Each number of the matrix is the synaptic weight from the neuron with label (column number -1) to neuron with label (row number -1).
