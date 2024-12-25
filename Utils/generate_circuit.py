
from qiskit import QuantumCircuit
from qiskit.circuit import Gate
import numpy as np

def qft(n):
    """创建QFT量子电路，n为量子比特的数量"""
    circuit = QuantumCircuit(n)
    gate_list = []

    # 创建QFT的逻辑
    for j in range(n):
        # 对每个量子比特进行Hadamard操作
        circuit.h(j)

        # 对每一对量子比特添加控制相位门
        for k in range(j + 1, n):
            theta = np.pi / float(2 ** (k - j))  # 相位角
            circuit.cp(theta, k, j)
            gate_list.append([k, j])  # 记录双量子门

    # 量子比特的交换操作
    for i in range(n // 2):
        circuit.swap(i, n - i - 1)
        gate_list.append([i, n - i - 1])  # 记录交换门

    return circuit, gate_list


def qft_with_cnot(n):
    """创建QFT电路，并将控制相位门转换为CNOT和Rz门"""
    circuit = QuantumCircuit(n)
    gate_list = []

    # QFT的实现，替换控制相位门为CNOT和Rz
    for j in range(n):
        # 对每个量子比特进行Hadamard操作
        circuit.h(j)

        # 替换控制相位门为CNOT和Rz
        for k in range(j + 1, n):
            theta = np.pi / float(2 ** (k - j))  # 相位角
            # 第一步：应用CNOT门
            circuit.cx(j, k)
            gate_list.append([j, k])  # 记录CNOT门

            # 第二步：应用Rz旋转门
            circuit.rz(theta, k)

    # 量子比特的交换操作
    for i in range(n // 2):
        circuit.swap(i, n - i - 1)
        gate_list.append([i, n - i - 1])  # 记录交换门

    return circuit, gate_list


if __name__ == '__main__':
    # 测试代码
    n = 100  # 比特数
    circuit, gate_list = qft_with_cnot(n)

    # 输出QFT电路和双量子门列表
    # print("Quantum Circuit:")
    # print(circuit)

    print("\nGate List (double-qubit gates):")
    print(gate_list)

    print(len(gate_list))

