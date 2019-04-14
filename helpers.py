from pyquil.noise import add_decoherence_noise
from pyquil.gates import *
from pyquil.quil import Program
import random
random.seed()
import numpy as np

pi = np.pi


def get_one_q_circuit(q_index, depth):
    """
    :param q_index: index of the qubit which the circuit acts on
    :depth: depth of the circuit
    :return: a program corresponding to a random U 
    """
    gate_set = [RX, RZ, T]
    instructions = []
    for i in range(depth):
        g = random.choice(gate_set)
        if g is T:
            instructions.append(RZ(pi/4,q_index))
        else:
            instructions.append(g(pi/2,q_index))
    
    return Program(instructions)

def get_two_q_circuit(q_index,n_cycles):
    """
    :param q_index: indexes of the qubits which the circuit acts on
    :n_cycles: depth of the circuit
    :return: a program corresponding to a random U 
    """
    get_set = [RX, RZ, T]
    instructions = []
    #1. applying Hadamard's in native language
    instructions.extend([RZ(pi/2, q_index[0]),RX(pi/2, q_index[0]),RZ(pi/2, q_index[0])])
    instructions.extend([RZ(pi/2, q_index[1]),RX(pi/2, q_index[1]),RZ(pi/2, q_index[1])])
    #2. applying CZ followed by 1 qubit gates 
    for i in range(n_cycles):
        instructions.append(CZ(q_index[0],q_index[1]))
        for idx in (q_index):
            g = random.choice(get_set)
            if g is T:
                instructions.append(RZ(pi/4,idx))
            else:
                instructions.append(g(pi/2,idx))
    
    return Program(instructions)


def add_pragma_block(program):
    inst = program.instructions
    new_inst = ['PRAGMA PRESERVE_BLOCK'] + inst + ['PRAGMA END_PRESERVE_BLOCK']
    return Program(new_inst)
    
def get_zx_DD_sequence(q_index, n):
    """
    :param q_index: index(es) of qubit(s) for applying DD sequence
    :param n:  number of sequence; each sequence is consisted of ZXZX pulses
    :return: program with DD sequence 
    """
    
    indexes = q_index
    if type(q_index) == int:
        q_index = [q_index]
    dd = []  
    for i, index in enumerate(q_index):
        dd.extend([RZ(pi, index),RX(pi,index), RZ(pi,index),RX(pi,index)] * n) #it can be modified to include buffer time (I gates)
            
        
    return Program(dd)


def get_xy_DD_sequence(q_index, n):
    """
    :param q_index: index(es) of qubit(s) for applying DD sequence
    :param n:  number of sequence; each sequence is consisted of XYXY (XY== RX(pi)RZ(pi)RX(pi)) pulses
    :return: program with DD sequence 
    """
    
    indexes = q_index
    if type(q_index) == int:
        q_index = [q_index]
    dd = []  
    for i, index in enumerate(q_index):
        dd.extend([RX(pi,index),RZ(pi, index),RX(pi,index),RX(pi,index), RZ(pi,index),RX(pi,index)] * n) 
            
        
    return Program(dd)
    
def get_idle_sequence(q_index, n, nI = 4):
    """
    :param q_index: index(es) of qubit(s) for applying DD sequence
    :param n:  number of wait circuits; each circuit consists of nI identity gates  
    :param nI: number of identity gates in wait circuit
    :return: program with wait sequence 
    """
    
    indexes = q_index
    if type(q_index) == int:
        q_index = [q_index]
    dd = []  
    for i, index in enumerate(q_index):
        dd.extend([I(index)] * (n * nI)) 
            
        
    return Program(dd)


# sampling programs with different gate times
def run_with_gate_time_sampling(cxn: QVMConnection,
                                programs: Iterable[Tuple[float, Program]],
                                program_modifier=None,
                                trials=20000):
    records = []
    base = 50e-9
    gate_times = np.array([1, 2, 3, 4]) * base

    for param, program in programs:
        program = program.copy()
        ro = program.declare('ro', 'BIT', 2)
        for gate_time in gate_times:
            noisy = add_decoherence_noise(program, gate_time_1q=gate_time, gate_time_2q=3 * gate_time).inst([
                MEASURE(0, ro[0]),
                MEASURE(1, ro[1]),
            ])

            if program_modifier:
                noisy = program_modifier(noisy)

            bitstring = np.array(cxn.run(noisy, [0, 1], trials))
            z0, z1 = np.mean(bitstring, axis=0)
            zz = 1 - (np.sum(bitstring, axis=1) % 2).mean() * 2

            f0, f1 = (trials - np.sum(bitstring, axis=0)) / trials
            ff = np.sum(np.sum(bitstring, axis=1) == 0) / trials

            record = {
                'z0': z0,
                'z1': z1,
                'zz': zz,
                'f0': f0,
                'f1': f1,
                'ff': ff,
                'param': param,
                'noise_param': gate_time,
            }
            records += [record]

    return records


# Computing mittigated values
def get_analyzed_and_mitigated(records):
    df_all = pd.DataFrame(records)

    noise_params = df_all['noise_param'].unique()

    qubits = 2
    mitigated = []

    for order in range(2, len(noise_params) + 1):
        matrix = noise_params[:order, np.newaxis] ** np.arange(order)

        mo = [[] for _ in range(qubits+1)]

        for param in df_all['param'].unique():
            df = df_all.query('{} == @{}'.format('param', 'param'))

            q1 = np.linalg.solve(matrix, df['z0'][:order])
            q2 = np.linalg.solve(matrix, df['z1'][:order])

            ff = np.linalg.solve(matrix, df['ff'][:order])

            mo[0] += [q1[0]] * len(df)
            mo[1] += [q2[0]] * len(df)
            mo[2] += [ff[0]] * len(df)

        mitigated += [mo]

    for order, o_values in enumerate(mitigated):
        for qubit, q_values in enumerate(o_values[:-1]):
            df_all.loc[:, 'm{}-{}'.format(qubit+1, order+1)] = np.array(q_values)
            df_all.loc[:, 'mf{}-{}'.format(qubit+1, order+1)] = 1 - np.array(q_values)
        df_all.loc[:, 'mfzz-{}'.format(order+1)] = np.array(o_values[-1])

    return df_all


# appling DD to a program
def add_dd(program: Program):
    new_program = program.copy_everything_except_instructions()

    counts = [0, 0]
    for gate in program:
        try:
            if len(gate.qubits) > 1:
                if abs(counts[0] - counts[1]) >= 2:
                    min_ind = int(counts[0] > counts[1])
                    times = max(int(abs(counts[0] - counts[1])/4), 1)

                    p = add_decoherence_noise(Program(get_dd_sec(min_ind)*times))

                    new_program.inst(p)
                counts = [0, 0]
            else:
                counts[gate.qubits[0].index] += 1
        except AttributeError:
            pass

        new_program.inst(gate)
    return new_program


# Generate Random cirquit
def two_qubit_circuit(length: int, qubit_one: int, qubit_two: int):
    """
    genereates two qubit identity equal circuit with given length

    :param length: length of the circuit
    :param qubit_one: one of the qubits
    :param qubit_two: second qubit
    :return: pyquil Program
    """

    p = Program()

    for j in range(int(length/2)):
        theta = 2 * np.pi * random.random()
        gate_list = [RZ(theta, qubit_one), RX(np.pi / 2, qubit_one), RX(- np.pi / 2, qubit_one),
                     CZ(qubit_one, qubit_two),
                     RZ(theta, qubit_two), RX(np.pi / 2, qubit_two), RX(- np.pi / 2, qubit_two), CZ(qubit_two, qubit_one)]
        new_gate = random.choice(gate_list)
        p.inst(new_gate)

    p += p.dagger()

    return Program('PRAGMA PRESERVE_BLOCK') + p + Program('PRAGMA END_PRESERVE_BLOCK')
