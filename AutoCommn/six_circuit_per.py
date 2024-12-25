import random

from qiskit import QuantumCircuit
import numpy as np
from Test.min_global_gate_num_test import statistics_gate_labels, count_gg_num
from Utils.division_by_irregularity import direct_calculation_of_tc_look_ahead, \
    transfer_qubit_list_by_gate_list_based_look_ahead, transfer_qubit_list_by_gate_list_based_look_ahead2, \
    split_into_k_parts
from Utils.generate_circuit import qft_with_cnot, qft
from Utils.transmission_cost_calculation_more import count_transfer_queue, is_global_gate
from qiskit.circuit import Gate, Parameter

import sys
sys.setrecursionlimit(50000)  # 设置更高的递归深度限制


def direct_calculation_of_tc_look_ahead_six(gate_list,cut_list):

    # 全局门标签
    statistics_gate_labels_list = statistics_gate_labels(gate_list, cut_list)


    # 量子位传输列表 ['q0', 'q0', 'q2', 'q0', 'q4', 'q1', 'q0', 'q4']
    initial_transfer_qubit_list0 = transfer_qubit_list_by_gate_list_based_look_ahead(gate_list, statistics_gate_labels_list)
    initial_transfer_qubit_list1 = transfer_qubit_list_by_gate_list_based_look_ahead2(gate_list,
                                                                                    statistics_gate_labels_list)
    # 统计合并传输队列
    #print('######################################################################################################')
    initial_transfer_queue0 = count_transfer_queue(gate_list, cut_list, initial_transfer_qubit_list0,
                                                  statistics_gate_labels_list)
    initial_transfer_queue1 = count_transfer_queue(gate_list, cut_list, initial_transfer_qubit_list1,
                                                   statistics_gate_labels_list)
    st0 = len(initial_transfer_queue0) * 2
    st1 = len(initial_transfer_queue1) * 2

    print('初始传输代价：' + str(min(st0,st1)),end=' ')
    # print('######################################################################################################')

    if st0 < st1:
        return initial_transfer_queue0
    else:
        return initial_transfer_queue1


def count_sublist_greater_than_p(two_d_list, p):
    # 统计长度大于p的子列表数量
    count_greater_than_p = sum(1 for sublist in two_d_list if len(sublist) >= p)

    # 计算总的子列表数量
    total_sublists = len(two_d_list)

    # 计算占比
    if total_sublists == 0:  # 防止除以0的错误
        return 0.0

    return count_greater_than_p / total_sublists


def qaoa_with_cnot(n, p):
    """
    生成 QAOA 电路，n 为量子比特数量，p 为 QAOA 的层数（迭代次数）。
    该电路中，所有多比特量子门都为 CNOT 门。
    返回：电路和量子门列表，gate_list 记录 CNOT 双量子门。
    """
    circuit = QuantumCircuit(n)
    gate_list = []

    # 初始化量子比特为叠加态（应用 Hadamard 门）
    for qubit in range(n):
        circuit.h(qubit)

    # 迭代 p 次，每一轮包括 Problem Hamiltonian 和 Mixing Hamiltonian
    for _ in range(p):
        # Problem Hamiltonian: 使用 CNOT 门模拟相位扩展操作
        for i in range(n - 1):
            # 使用 CNOT 门模拟控制相位
            circuit.cx(i, i + 1)
            gate_list.append([i, i + 1])  # 记录 CNOT 双量子门

        # Mixing Hamiltonian: 继续使用 CNOT 门
        for i in range(n - 1):
            # 这里我们也使用 CNOT 门来混合态
            circuit.cx(i, i + 1)
            gate_list.append([i, i + 1])  # 记录 CNOT 双量子门

    return circuit, gate_list


def bv_algorithm(n, a):
    """
    生成BV（Bernstein-Vazirani）算法电路
    :param n: 量子比特数（即目标字符串a的长度）
    :param a: 二进制字符串，表示未知的常数向量a（例如 '101'）
    :return: 电路和量子门列表（包含CNOT门）
    """
    # 创建量子电路
    circuit = QuantumCircuit(n + 1, n)  # n + 1 个量子比特（额外的辅助比特），n个经典比特
    gate_list = []

    # 1. 初始化量子比特，先应用Hadamard门
    for qubit in range(n):
        circuit.h(qubit)  # 将n个比特都初始化为叠加态

    circuit.x(n)  # 将额外的辅助比特设置为 |1⟩
    circuit.h(n)  # 对辅助比特应用Hadamard变换

    # 2. 应用Oracle：模拟函数 f(x) = a · x mod 2
    # oracle 操作实现了 CNOT 和相应的控制相位操作
    for qubit in range(n):
        if a[qubit] == '1':  # 根据 a 的每个比特，决定是否应用 CNOT 门
            circuit.cx(qubit, n)  # 如果 a 的对应位为 1，则在量子比特 qubit 和辅助比特 n 之间应用 CNOT 门
            gate_list.append([qubit, n])  # 记录CNOT门

    # 3. 对所有输入量子比特再次应用 Hadamard 门
    for qubit in range(n):
        circuit.h(qubit)

    # 4. 测量输出量子比特（将结果存储在经典比特中）
    circuit.measure(range(n), range(n))

    # 打印电路
    # print(circuit.draw())

    # 检查 gate_list 中的门
    for gate in gate_list:
        # 检查是否是CNOT门
        if isinstance(circuit.data[gate[0]][0], Gate) and circuit.data[gate[0]][0].name == 'cx':
            print(f"CNOT Gate found at qubits {gate[0]} and {gate[1]}")

    return circuit, gate_list

def generate_random_a(length):
    # 随机生成一个二进制字符串
    a = ''.join(random.choice('01') for _ in range(length))
    print(a)
    return a


def uccsd_ansatz(n):
    """
    生成 UCCSD（Unitary Coupled Cluster）ansatz的量子电路，
    并将其转换为仅包含 CNOT 门的多量子比特电路。

    :param n: 量子比特数
    :return: 量子电路和 CNOT 门列表
    """
    circuit = QuantumCircuit(n)
    gate_list = []

    # 生成 singles 操作：选择所有不同的 (i, j) 进行单激发
    singles = [(i, j) for i in range(n) for j in range(i + 1, n)]

    # 生成 doubles 操作：选择所有不同的 (i, j, k, l) 进行双激发
    doubles = [(i, j, k, l) for i in range(n) for j in range(i + 1, n) for k in range(j + 1, n) for l in
               range(k + 1, n)]

    # 1. 初始化量子比特，使用Hadamard门
    for qubit in range(n):
        circuit.h(qubit)

    # 2. 添加 singles 激发操作：CNOT 门
    for (i, j) in singles:
        # 在 i 和 j 之间添加 CNOT 门
        circuit.cx(i, j)
        gate_list.append((i, j))  # 记录 CNOT 门

    # 3. 添加 doubles 激发操作：模拟双激发操作
    for (i, j, k, l) in doubles:
        # 在 (i, j) 和 (k, l) 之间添加 CNOT 门（双激发）
        circuit.cx(i, j)
        circuit.cx(k, l)
        gate_list.append((i, j))
        gate_list.append((k, l))

    return circuit, gate_list


def toffoli_to_cnot(circuit, control1, control2, target, gate_list):
    """
    将 Toffoli 门转换为 5 个 CNOT 门，并将其加入到 gate_list 中
    :param circuit: 量子电路
    :param control1: 控制比特 1
    :param control2: 控制比特 2
    :param target: 目标比特
    :param gate_list: 用于存储 CNOT 门的列表
    :return: 更新后的量子电路和门列表
    """
    # 第一步：应用 CNOT 门，控制比特1 -> 目标比特
    circuit.cx(control1, target)
    gate_list.append([control1, target])

    # 第二步：应用 CNOT 门，控制比特2 -> 目标比特
    circuit.cx(control2, target)
    gate_list.append([control2, target])

    # 第三步：应用 CNOT 门，控制比特1 -> 控制比特2
    circuit.cx(control1, control2)
    gate_list.append([control1, control2])

    # 第四步：应用 CNOT 门，控制比特2 -> 目标比特
    circuit.cx(control2, target)
    gate_list.append([control2, target])

    # 第五步：应用 CNOT 门，控制比特1 -> 目标比特
    circuit.cx(control1, target)
    gate_list.append([control1, target])

    return circuit, gate_list


def ripple_carry_adder(n):
    """
    构建一个 Ripple-Carry Adder (RCA) 量子电路
    :param n: 量子比特数
    :return: 量子电路和 CNOT 门列表
    """
    # 创建量子电路，包含两个输入数A和B和一个进位比特
    circuit = QuantumCircuit(n * 2 + 1, n)
    gate_list = []

    # 初始化输入量子比特 A 和 B（假设它们已经是给定的二进制数）
    # 假设 A 和 B 的比特已经在量子比特中，通过一些初始门实现

    # 计算逐位加法，生成进位并处理
    for i in range(n):
        # 计算 A[i] + B[i] + carry
        if i == 0:
            circuit.cx(i, n + i)  # A[i] 和 B[i] 进行加法
            gate_list.append([i, n + i])  # 记录 CNOT 门

        else:
            # 逐位加法，处理进位
            toffoli_to_cnot(circuit, i, n + i, n + i + 1, gate_list)

    return circuit, gate_list


from qiskit import QuantumCircuit


def mct_decompose(circuit, n):
    """
    递归分解多控制非门 (MCT gate) 到 Toffoli 门和 CNOT 门。

    参数：
        circuit (QuantumCircuit): 量子电路对象
        n (int): 控制比特的数量（包括目标比特在内的量子比特数量）
    """
    # 基本情况：只有 2 个控制比特时，直接用 Toffoli 门
    if n == 2:
        controls = circuit.qregs[0][:2]  # 取前两个比特作为控制比特
        target = circuit.qregs[0][2]  # 目标比特
        circuit.toffoli(controls[0], controls[1], target)
        return

    # 递归情况：n > 2 时，递归分解 MCT
    controls = circuit.qregs[0][:n - 1]  # 获取前 n-1 个控制比特
    target = circuit.qregs[0][n - 1]  # 获取最后一个比特作为目标比特

    # 创建辅助量子比特，用于递归拆解
    ancillas = [circuit.qregs[0][i] for i in range(n - 2)]  # 辅助比特的数量是 n-2

    # 递归分解 MCT 门
    # 步骤 1：将 MCT 门分解成若干个 Toffoli 和 CNOT 门
    for i in range(n - 2):
        # 使用 CNOT 门将控制比特和辅助比特连接
        circuit.cx(controls[i], ancillas[i])

    # 步骤 2：使用 Toffoli 门进行递归分解
    mct_decompose(circuit, n - 1)

    # 步骤 3：再次连接 CNOT 门
    for i in range(n - 2):
        circuit.cx(controls[i], ancillas[i])

    # 步骤 4：最终使用 Toffoli 门完成目标操作
    circuit.toffoli(controls[0], controls[1], target)


if __name__ == '__main__':

    # circuit, gate_list = qft(100)

    # circuit, gate_list = qaoa_with_cnot(300, 20)

    # num_list = [99,199,299,399,499,599]

    node_list = [2,10,20,50,100]

    result = []

    for node in node_list:
        num = 299
        circuit, gate_list = bv_algorithm(num, generate_random_a(num))

        # circuit, gate_list = uccsd_ansatz(16)

        # circuit, gate_list = ripple_carry_adder(150)
        # circuit = QuantumCircuit(10)
        # mct_decompose(circuit, 10)
        # print(circuit)

        print(len(gate_list))

        cut_list = split_into_k_parts(num + 1, node)

        global_count = is_global_gate(gate_list, cut_list)[0].count(1) * 2
        print(global_count)

        transfer_queue = direct_calculation_of_tc_look_ahead_six(gate_list, cut_list)

        tc = len(transfer_queue) * 2
        print(tc)

        print(global_count / tc)
        result.append(global_count / tc)
        # print(transfer_queue)
        # for i in range(30):
        #     print(count_sublist_greater_than_p(transfer_queue, i+1))

    for i in result:
        print(i)