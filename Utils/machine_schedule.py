'''
功能： 根据待分布式的线路和ibm提供的量子计算机，选择可用的量子计算机，输出所有满足运行条件的量子计算机选择
ibm共提供21台量子计算机，量子位数分别为：127、65、27*9、16、7*5、5*4
本代码主要面向于线路划分为两部分
'''
import re
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import Utils.toffoli as utils
import itertools

'''读取qasm文件并进行存储'''
def converter_circ_from_qasm(input_file_name):
    gate_list = []
    qbit = 0  # 量子位
    qasm_file = open(input_file_name, 'r')
    iter_f = iter(qasm_file)
    reserve_line = 0
    num_line = 0
    for line in iter_f:  # 遍历文件，一行行遍历，读取文本
        num_line += 1
        if num_line <= reserve_line:
            continue
        else:
            if 'qreg' in line:
                qbit = get_data(line)[0]
            if line[0:1] == 'x' or line[0:1] == 'X':
                '''获取X门'''
                x = get_data(line)
                x_target = x[0]
                listSingle = [x_target]
                gate_list.append(listSingle)
            if line[0:2] == 'CX' or line[0:2] == 'cx':
                '''获取CNOT'''
                cnot = get_data(line)
                cnot_control = cnot[0]
                cnot_target = cnot[1]
                listSingle = [cnot_control, cnot_target]
                gate_list.append(listSingle)
            if line[0:2] == 'CP' or line[0:2] == 'cp':
                cp = get_data(line)
                cp_one = cp[1]
                cp_two = cp[2]
                listSingle = [cp_one,cp_two]
                gate_list.append(listSingle)
            if line[0:4] == 'SWAP' or line[0:4] == 'swap':
                swap = get_data(line)
                swap_one = swap[0]
                swap_two = swap[1]
                cnot_one = [swap_one,swap_two]
                cnot_two = [swap_two,swap_one]
                gate_list.append(cnot_one)
                gate_list.append(cnot_two)
                gate_list.append(cnot_one)
            if line[0:3] == 'CCX' or line[0:3] == 'ccx':
                '''获取toffoli'''
                toffoli = get_data(line)
                toffoli_control1 = toffoli[0]
                toffoli_control2 = toffoli[1]
                toffoli_target = toffoli[2]
                listSingle = [toffoli_control1, toffoli_control2, toffoli_target]
                gate_list.append(listSingle)
    return gate_list, qbit

def get_data(str):
    pattern = re.compile("[\d]+")
    result = re.findall(pattern, str)
    return result

'''
将gate_list全部转换为int
'''
def list_str_to_int(gate_list):
    new_gate_list = []
    for i in range(len(gate_list)):
        son_new_gate_list = list(map(int, gate_list[i]))
        new_gate_list.append(son_new_gate_list)
    return new_gate_list

'''gate_list去除单门'''
def remove_single_qubit_gate(gate_list):
    i = 0
    while i<len(gate_list):
        if len(gate_list[i]) == 1 :
            gate_list.pop(i)
            i -= 1
        i += 1
    return gate_list

'''输出所有可能的物理架构  异构'''
def all_results(qubit_list,ciucuit_qbit):
    j = 0
    first_part = []
    second_part = []
    result_list = []
    for i in range(len(qubit_list)):
        for j in range(i, len(qubit_list)):
            if i == j:
                continue
            first_part = qubit_list[i]
            second_part = qubit_list[j]
            if first_part + second_part >= circuit_qbit + 1:  #a1+a2>=qubit+1
                if first_part < circuit_qbit and second_part < circuit_qbit:  #a1<qubit a2<qubit
                    result_list.append([first_part, second_part])
                else:
                    continue
    # print(result_list)
    new_result_list = []
    for k in result_list:
        if k not in new_result_list:
            new_result_list.append(k)  #去掉列表里重复的
     #new_result_list:所有可能的物理架构
    return new_result_list


'''找出同构'''
def same_architecture(new_result_list):
    same_architecture_list = []
    for i in range(len(new_result_list)):
        if new_result_list[i][0] == new_result_list[i][1]: #a1=a2同构
            same_architecture_list.append(new_result_list[i])
    return same_architecture_list


'''本源架构'''
def by_architecture(by_qubit,circuit_qbit):
    byuan_architecture_list = all_results(by_qubit,circuit_qbit)
    return byuan_architecture_list


'''门列表转成邻接矩阵'''
def matrixtable(gate_list, circuit_qbit):
    matrix = np.zeros((circuit_qbit, circuit_qbit))
    for gate in gate_list:
        for j in range(len(gate)-1):
            matrix[gate[j]][gate[j+1]] = 1+matrix[gate[j]][gate[j+1]]
    return matrix

# 绘制图
def adMmatrix2Img(matrix):
    G = nx.Graph()
    n = len(matrix)
    point = []
    for i in range(n):
        point.append(i)
    G.add_nodes_from(point)
    edglist = []
    for i in range(n):
        for k in range(i + 1, n):
            if matrix[i][k] > 0:
                edglist.append((i, k))
                G.add_edge(i,k,weight = int(matrix[i][k]+matrix[k][i]))
    # G.add_edges_from(edglist)
    position = nx.circular_layout(G)
    weights = nx.get_edge_attributes(G, "weight")
    nx.draw_networkx_nodes(G, position, nodelist=point, node_color="r")
    nx.draw_networkx_edges(G, position)
    nx.draw_networkx_labels(G, position)
    nx.draw_networkx_edge_labels(G, position, edge_labels = weights)
    plt.show()


'''
输出所有种线序组合
'''
def line_sequence_change_combination(qbit,cut_point):
    str1 = ''
    # 最多支持 26量子位
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z']
    for i in range(int(qbit)):
        # str1：ABCDEFGHIJKL.... 共qbit个
        str1 = str1 + str2[i]
    # 排列组合 将str1按qibt/2 一分为，共有C(qbit/2,qbit)种情况
    line_sequence_combination = []  # 线序排列集合['ABCDEF', 'ABCDEG', 'ABCDEH', 'ABCDEI', 'ABCDEJ', 'ABCDEK', 'ABCDEL', 'ABCDFG', 'ABCDFH', 'ABCDFI', 'ABCDFJ'......]
    for i in itertools.combinations(str1, cut_point):
        #  print(''.join(i), end=" ")
        line_sequence_combination.append(''.join(i))
    return line_sequence_combination

'''计算全局门个数 （门列表，切线位置，量子位数）'''
def count_global_gate_num(gate_list,cut_point):
    global_gate_list = []
    local_gate_list=[]
    count_num = 0
    cut_list_one = []
    cut_list_two = []
    # 第一个分区的所有量子位集合
    for i in range(cut_point):
        cut_list_one.append(i)
    # print(cut_list_one)
    # 第二个分区的所有量子位集合
    for i in range(cut_point,circuit_qbit):
        cut_list_two.append(i)
    # print(cut_list_two)
    # 门列表与分区的交集为门列表本身，并且与另一个分区没有交集，是局部门，否则为全局门
    for i in range(len(gate_list)):
        res1 = list(set(gate_list[i]) & set(cut_list_one))
        res2 = list(set(gate_list[i]) & set(cut_list_two))
        if (len(res1)==len(set(gate_list[i])) and (len(res2)==0)) or (len(res2)==len(set(gate_list[i])) and (len(res1)==0)):
            local_gate_list.append(gate_list[i])
        else:
            global_gate_list.append(gate_list[i])
            count_num = count_num + 1
    # print(global_gate_list)
    # print(local_gate_list)
    return count_num

'''
根据前半部分线序推出全部线序  ABCDFG =>ABCDFGEHIJKL
'''
def single_to_all_line_sequence(line_sequence_combination,qbit):
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
            'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
    all_line_sequence = list(line_sequence_combination)
    for i in range(int(qbit)):
        if str2[i] not in all_line_sequence:
            all_line_sequence.append(str2[i])
    new_all_line_sequence = "".join(all_line_sequence)
    return new_all_line_sequence

'''根据线序改变gate_list'''
def change_gate_list_by_line_sequence(initial_line_sequence,gate_list,single_line_sequence):
    new_gate_list = list_str_to_int(gate_list)
    all_line_sequence = single_to_all_line_sequence(single_line_sequence, circuit_qbit)  # ABCDEG FHIJKL
    print(all_line_sequence)
    for j in range(len(all_line_sequence)):  # j: 0-11
        if all_line_sequence[j] == initial_line_sequence[j]:
            continue
        if all_line_sequence[j] != initial_line_sequence[j]:
            for k in range(len(new_gate_list)):  # k: 门数
                for p in range(len(new_gate_list[k])):  # l:0-3
                    if new_gate_list[k][p] == letter_to_number(all_line_sequence[j]):
                        new_gate_list[k][p] = str(j)
    new_gate_list = list_str_to_int(new_gate_list)
    return new_gate_list

'''
把字母变成对于的数字
'''
def letter_to_number(letter):
    number = ord(letter) - 65
    return number

'''判断一个门是否为全局门'''
def judge_is_global_gate(gate,cut_point):
    cut_list_one = []
    cut_list_two = []
    # 第一个分区的所有量子位集合
    for i in range(cut_point):
        cut_list_one.append(i)
    # 第二个分区的所有量子位集合
    for i in range(cut_point, circuit_qbit):
        cut_list_two.append(i)
    res1 = list(set(gate) & set(cut_list_one))
    res2 = list(set(gate) & set(cut_list_two))
    if (len(res1) == len(set(gate)) and (len(res2) == 0)) or (
            len(res2) == len(set(gate)) and (len(res1) == 0)):
        is_global_gate = 0
    else:
        is_global_gate = 1
    return is_global_gate


'''
根据量子位数、割点、门列表，最后计算各划分方法下的最低量子门数或传输代价
'''
def count_transfer_queue(cut_point,gate_list):
    transfer_queue = [] # 传输队列 二维数组 每个子数组表示一次合并传输 传输队列长度*2即为传输代价
    # 生成全局门列表
    global_gate_list = []
    for gate in gate_list:
        # 是全局门
        if judge_is_global_gate(gate,cut_point) == 1:
            global_gate_list.append(gate)
    # 循环遍历每个全局门 寻找此全局门往后的全局门能否和此门合并传输
    while len(global_gate_list) > 0 :
        megred_trasfer_list = []  # 合并传输列表
        current_global_gate = global_gate_list[0] # 当前全局门
        megred_trasfer_list.append(global_gate_list[0])
        j = 1
        while 0 < len(global_gate_list):
            if j < len(global_gate_list):
                next_global_gate = global_gate_list[j]  # 下一个全局门
            # 判断是否满足合并传输  具有相同的量子位且两门之间无间隔门或间隔门不影响
            if len(megred_trasfer_list) == 1: # 只有一个门
                same_qubit = list(set(current_global_gate) & set(next_global_gate))  # 公共量子位
                if same_qubit != []:  # 集合交集非空，则存在相同元素，具有公共量子位 且传输方向一致
                    if gate_list.index(current_global_gate) + 1 == gate_list.index(next_global_gate) :
                        megred_trasfer_list.append(next_global_gate)  # 将此门加入合并传输列表中
                        j += 1
                        continue
                else: #无法合并传输
                    break
            if len(megred_trasfer_list) > 1: # 两个门以上
                same_qubit = judge_megred_transfer(megred_trasfer_list) # 公共量子位
                # 当前门与合并传输列表中的门具有相同量子位
                if list(set(next_global_gate) & set(same_qubit)) != []:
                    # 这些门都相邻
                    if gate_list.index(next_global_gate) == gate_list.index(megred_trasfer_list[-1]) +1 :
                        megred_trasfer_list.append(next_global_gate)  # 将此门加入合并传输列表中
                        j += 1
                        continue
                else: #无法合并传输
                    break
            # 循环停止标志
            if j >= len(global_gate_list):
                break
        transfer_queue.append(megred_trasfer_list) # 传输队列中加入合并传输列表
        # 从global_gate_list中删除已经合并传输的门
        for j in range(len(megred_trasfer_list)):
            global_gate_list.pop(0)
        print(transfer_queue)
    return transfer_queue


'''判断两门是否能够合并传输，并返回合并传输的量子位'''
def judge_megred_transfer(megred_trasfer_list):
    same_qubit = list(set(megred_trasfer_list[0]) & set(megred_trasfer_list[1]))  # 前两个门公共量子位
    for i in range(2,len(megred_trasfer_list)):
        same_qubit = list(set(same_qubit) & set(megred_trasfer_list[i]))
    return same_qubit


'''判断这个门的公共量子位所在分区存在几个量子位'''
def count_qubit_in_one_partation(same_qubit,gate,cut_point):
    count_qubit_num = 0
    cut_list_one = []
    cut_list_two = []
    # 第一个分区的所有量子位集合
    for i in range(cut_point):
        cut_list_one.append(i)
    # 第二个分区的所有量子位集合
    for i in range(cut_point, circuit_qbit):
        cut_list_two.append(i)
    # print(cut_list_one)
    # print(cut_list_two)
    # 判断分区是否仅存在一个量子位
    if same_qubit in cut_list_one:
         count_qubit_num = len(list(set(cut_list_one)&set(gate)))
    if same_qubit in cut_list_two:
         count_qubit_num = len(list(set(cut_list_two)&set(gate)))
    if len(same_qubit) == 3:
        count_qubit_num = 1
    return count_qubit_num

'''判断传输方向是否一致1.0'''
def judge_transfer_direction(current_global_gate,next_global_gate):
    same_qubit = list(set(current_global_gate) & set(next_global_gate))
    is_same_direction = 0 # 0表示传输方向不一致
    if len(current_global_gate) == 3: #toffoli
        # 公共量子位所在分区存在两个量子位
        if count_qubit_in_one_partation(same_qubit,current_global_gate,cut_point) == 2:
            direction_of_current_global_gate = 1 # 反向传输
            is_same_direction = 0  # 0表示传输方向不一致
        # 公共量子位所在分区存在1个量子位
        if count_qubit_in_one_partation(same_qubit, current_global_gate, cut_point) == 1:
            direction_of_current_global_gate = 0  # 正向向传输
            if len(next_global_gate) == 3: #toffoli
                if count_qubit_in_one_partation(same_qubit,next_global_gate,cut_point) == 2:
                    is_same_direction = 0  # 0表示传输方向不一致
                if count_qubit_in_one_partation(same_qubit,next_global_gate,cut_point) == 1:
                    is_same_direction = 1  # 1表示传输方向一致
            if len(next_global_gate) == 2:  # cnot
                is_same_direction = 1
    if len(current_global_gate) == 2:  # cnot
        if len(next_global_gate) == 3:  # toffoli
            if count_qubit_in_one_partation(same_qubit, next_global_gate, cut_point) == 2:
                is_same_direction = 0  # 0表示传输方向不一致
            if count_qubit_in_one_partation(same_qubit, next_global_gate, cut_point) == 1:
                is_same_direction = 1  # 1表示传输方向一致
        if len(next_global_gate) == 2:  # cnot
            is_same_direction = 1  # 1表示传输方向一致
    return is_same_direction



'''判断传输方向是否一致2.0'''
def judge_transfer_direction_2(current_global_gate,next_global_gate):
    is_same_direction = 0  # 0表示传输方向不一致
    same_qubit = list(set(current_global_gate) & set(next_global_gate))
    if count_qubit_in_one_partation(same_qubit, current_global_gate, cut_point) == 1 and count_qubit_in_one_partation(same_qubit, next_global_gate, cut_point) == 1:
        is_same_direction = 1
    return is_same_direction

if __name__ == '__main__':
    # 读取线路
    input_filename = 'E:/python/Distributed_quantum_circuits_scheduling/qasm/ham7_104.qasm'
    # gate_list = converter_circ_from_qasm(input_filename)[0]
    circuit_qbit = int(converter_circ_from_qasm(input_filename)[1])  # 量子位
    print('量子位数：' + str(circuit_qbit))
    gate_list = list_str_to_int(converter_circ_from_qasm(input_filename)[0])
    print('qasm读取线路:', end=' ')
    print(gate_list)
    gate_list = remove_single_qubit_gate(gate_list)
    print(gate_list)
    # decompose_gate_list = utils.toffoli_decompose_c(gate_list)
    # print('toffoli门分解后线路:', end=' ')
    # print(decompose_gate_list)

    # 选择架构，输出所有可选架构
    ibm_qubit = [127, 65, 27, 27, 16, 7, 7, 5, 5]  # ibm架构的量子位数
    byuan_qubit = [6, 6]                           # 本源量子架构的量子位数
    # 所有可能的物理架构 异构
    new_qubit_result_list = all_results(ibm_qubit,circuit_qbit)
    print("所有可能的物理架构（异构）:",end='')
    print(new_qubit_result_list)
    # 所有可能的物理架构 同构
    same_architecture_list = same_architecture(new_qubit_result_list)
    print("所有可能的物理架构（同构）:", end='')
    print(same_architecture_list)
    # 所有可能的物理架构 本源量子
    by_architecture_list = by_architecture(byuan_qubit,circuit_qbit)
    print("所有可能的物理架构（本源）:", end='')
    print(by_architecture_list)

    # # 门列表转成邻接矩阵
    # gate_matrixtable = matrixtable(decompose_gate_list, circuit_qbit)
    # print("邻接矩阵:")
    # print(gate_matrixtable)
    # # 根据矩阵绘制图像
    # adMmatrix2Img(gate_matrixtable)

    # cut_point = 5
    # # 计算全局门的个数
    # count_num = count_global_gate_num(gate_list,cut_point)
    # print("初始全局门的个数：", end='')
    # print(count_num)
    #
    # line_sequence_combination = line_sequence_change_combination(circuit_qbit,cut_point)
    # print(line_sequence_combination)
    # initial_line_sequence = single_to_all_line_sequence(line_sequence_combination[0], circuit_qbit)
    # # new_gate_list = change_gate_list_by_line_sequence(initial_line_sequence,gate_list,'ABD')
    # # print(count_global_gate_num(new_gate_list,3))
    # # 遍历每一种线序下
    # count_global_gate_num_list = []
    # for i in range(len(line_sequence_combination)):
    #     # print(line_sequence_combination[i])
    #     new_gate_list = change_gate_list_by_line_sequence(initial_line_sequence, gate_list, line_sequence_combination[i])
    #     print(new_gate_list)
    #     print(count_global_gate_num(new_gate_list,cut_point))
    #     count_global_gate_num_list.append(count_global_gate_num(new_gate_list,cut_point))
    #     print('')
    # print(count_global_gate_num_list)
    # print('最小全局门数：',end='')
    #
    # min_global_gate_num = min(count_global_gate_num_list)
    # print(min_global_gate_num)
    #
    # print('###################################################################')
    # # 输出最小全局门数下的 gate_list
    # for i in range(len(count_global_gate_num_list)):
    #     if count_global_gate_num_list[i] == min_global_gate_num:
    #         min_line_sequence = line_sequence_combination[i]
    #         min_global_gate_num_gate_list = change_gate_list_by_line_sequence(initial_line_sequence, gate_list, line_sequence_combination[i])
    #         print(min_global_gate_num_gate_list)
    #         print(count_global_gate_num(min_global_gate_num_gate_list,cut_point))
    # print('###################################################################')
    #
    # transfer_queue = count_transfer_queue(cut_point, gate_list)
    # print('###################################################################')
    # print(transfer_queue)
    # print(len(transfer_queue))
    # print('###################################################################')
    # cut_point = 2
    # print(judge_transfer_direction([0,1,2],[1,2]))
    # print(judge_transfer_direction_2([0, 1, 2], [1, 2]))

