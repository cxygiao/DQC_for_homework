from qiskit import Aer
from qiskit.aqua import QuantumInstance
from qiskit.aqua.algorithms import QAOA
from qiskit.aqua.components.optimizers import SPSA
from qiskit.optimization.applications.ising import graph_partition

# 构建图的最小割问题
num_vertices = 4
edges = [(0, 1, 1), (0, 2, 2), (1, 2, 3), (1, 3, 4), (2, 3, 5)]
qubit_op, offset = graph_partition.get_operator(edges)

# 构建QAOA实例
optimizer = SPSA(max_trials=100)
qaoa = QAOA(qubit_op, optimizer, p=1)

# 运行QAOA算法并获取结果
backend = Aer.get_backend('statevector_simulator')
quantum_instance = QuantumInstance(backend)
result = qaoa.run(quantum_instance)

# 输出结果
print("Minimum cut:", graph_partition.get_graph_solution(result.eigenstate))
print("Minimum cut value:", result.eigenvalue + offset)
