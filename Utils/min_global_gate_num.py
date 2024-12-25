'''
最小化全局门数量
陈新宇
2024.11.29
'''
import datetime
import re
import time
from multiprocessing.pool import Pool
import itertools
import random
import Utils.transmission_cost_calculation_more as TC
from Utils.Switching_line_sequences import list_str_to_int
from Utils.machine_schedule import single_to_all_line_sequence, letter_to_number

'''计算每个门的标签'''
def statistics_gate_labels(gate_list,cut_list):
    statistics_gate_labels_list = [] #门标签
    partition_interval_list = [] # 分区区间 [[0,1],[2,3],[4,6]]
    c = 0
    o = 0
    while c < len(cut_list):
        # 生成区间
        partition_interval_list.append([o,o+cut_list[c]-1])
        o = o+cut_list[c]
        c += 1
    for i in range(len(gate_list)): # 每个门
        gate_labels = []
        for j in range(len(gate_list[i])): # 每个门的每个量子位
            for k in range(len(partition_interval_list)):
                if gate_list[i][j] >= partition_interval_list[k][0] and gate_list[i][j] <= partition_interval_list[k][1]:
                    gate_labels.append(k)
        statistics_gate_labels_list.append(gate_labels)
    # print(statistics_gate_labels_list)
    return statistics_gate_labels_list

'''判断gate_list中门是不是全局门'''
def is_global_gate(gate_list,cut_list):
    is_global_gate_list = []
    global_gate_num_list = []
    statistics_gate_labels_list = statistics_gate_labels(gate_list,cut_list)
    for i in range(len(statistics_gate_labels_list)):
        # 局部门
        if statistics_gate_labels_list[i][0] == statistics_gate_labels_list[i][1]:
            is_global_gate_list.append(0)
        # 全局门
        else:
            is_global_gate_list.append(1)
            global_gate_num_list.append('g'+ str(i))
    return is_global_gate_list,global_gate_num_list

# 第二种方法计算全局门数量
'''判断一个门是否为全局门'''
def judge_is_global_gate(gate,cut_point):
    cut_list_one = []
    cut_list_two = []
    # 第一个分区的所有量子位集合
    for i in range(cut_point):
        cut_list_one.append(i)
    # 第二个分区的所有量子位集合
    for i in range(cut_point, circuit_qubit):
        cut_list_two.append(i)
    res1 = list(set(gate) & set(cut_list_one))
    res2 = list(set(gate) & set(cut_list_two))
    if (len(res1) == len(set(gate)) and (len(res2) == 0)) or (
            len(res2) == len(set(gate)) and (len(res1) == 0)):
        is_global_gate = 0
    else:
        is_global_gate = 1
    return is_global_gate

# 第二种方法计算全局门数量
def count_gg_num(gate_list,cut_point):
    count = 0
    for i in range(len(gate_list)):
        if judge_is_global_gate(gate_list[i],cut_point)==1:
            count+=1
    return count

'''根据线序计算全局门数'''
def count_gg_num_by_line_sequence(initial_line_sequence,gate_list,single_line_sequence):
    new_gate_list = change_gate_list_by_line_sequence(initial_line_sequence,gate_list,single_line_sequence)
    print(new_gate_list)
    gg_num = len(is_global_gate(new_gate_list,cut_list)[1])
    return gg_num

'''根据线序计算全局门数'''
def k_count_gg_num_by_line_sequence(initial_line_sequence,gate_list,line_sequence):
    new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,line_sequence)
    gg_num = len(is_global_gate(new_gate_list,cut_list)[1])
    return gg_num

'''根据线序改变gate_list(适用于两个分区)'''
def change_gate_list_by_line_sequence(initial_line_sequence,gate_list,single_line_sequence):
    qubit = max([max(row) for row in gate_list]) + 1
    new_gate_list = list_str_to_int(gate_list)
    all_line_sequence = single_to_all_line_sequence(single_line_sequence, qubit)  # ABCDEG FHIJKL
    # print(all_line_sequence)
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


'''根据线序改变gate_list(适用于K个分区)'''
def k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,all_line_sequence):
    new_gate_list = list_str_to_int(gate_list)
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


'''2分区下计算最小全局门数'''
def count_min_global_gate_num(gate_list,cut_list):
    qubit = max([max(row) for row in gate_list]) + 1
    str1 = ''
    # 最多支持 62量子位
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
            'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
    for i in range(int(qubit)):
        # str1：ABCDEFGHIJKL.... 共qbit个
        str1 = str1 + str2[i]
    initial_line_sequence = str1
    # 排列组合 将str1按qibt/2 一分为，共有C(qbit/2,qbit)种情况
    line_sequence_combination = []  # 线序排列集合['ABCDEF', 'ABCDEG', 'ABCDEH', 'ABCDEI', 'ABCDEJ', 'ABCDEK', 'ABCDEL', 'ABCDFG', 'ABCDFH', 'ABCDFI', 'ABCDFJ'......]
    min_gg_bum = len(is_global_gate(gate_list,cut_list)[1])
    min_line_sequence = initial_line_sequence[0:cut_list[0]]
    for i in itertools.combinations(str1, cut_list[0]):
        print(''.join(i), end=" ")
        # line_sequence_combination.append(''.join(i))
        line_sequence = (''.join(i))
        # all_line_sequence = single_to_all_line_sequence(line_sequence, qubit)  # ABCDEG FHIJKL
        gg_num = count_gg_num_by_line_sequence(initial_line_sequence,gate_list,line_sequence)
        print(gg_num)
        if gg_num < min_gg_bum:
            min_gg_bum = gg_num
            min_line_sequence = line_sequence
    min_gate_list = change_gate_list_by_line_sequence(initial_line_sequence,gate_list,min_line_sequence)
    return min_gg_bum,min_gate_list


'''K分区下计算最小全局门数'''
def k_count_min_global_gate_num(gate_list,cut_list,line_sequence_list):
    initial_line_sequence = ''
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
            'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    qubit = max([max(row) for row in gate_list]) + 1
    for i in range(int(qubit)):
        # str1：ABCDEFGHIJKL.... 共qbit个
        initial_line_sequence = initial_line_sequence + str2[i]
    min_gg_bum = len(is_global_gate(gate_list, cut_list)[1])
    min_line_sequence = initial_line_sequence
    for i in range(len(line_sequence_list)):
        line_sequence = line_sequence_list[i]
        gg_num = count_gg_num_by_line_sequence(initial_line_sequence, gate_list, line_sequence)
        if gg_num < min_gg_bum:
            min_gg_bum = gg_num
            min_line_sequence = line_sequence
    min_gate_list = change_gate_list_by_line_sequence(initial_line_sequence, gate_list, min_line_sequence)
    return min_gg_bum, min_gate_list


'''初始线路'''
def Initial_line_sequence(gate_list):
    initial_line_sequence = ''
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
            'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    qubit = max([max(row) for row in gate_list]) + 1
    for i in range(int(qubit)):
        # str1：ABCDEFGHIJKL.... 共qbit个
        initial_line_sequence = initial_line_sequence + str2[i]
    return initial_line_sequence


'''K分区下计算最小传输代价'''
def k_count_min_st_num2(gate_list,cut_list,line_sequence_list):
    initial_line_sequence = Initial_line_sequence(gate_list)
    ST = []
    STO = []
    lenl = len(line_sequence_list)
    print(lenl)
    minio = len(gate_list)*2
    mini = len(gate_list)*2
    for i in range(lenl):
        print(i, end=':')
        print(line_sequence_list[i],end='')
        new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,line_sequence_list[i])
        # print(new_gate_list,end=' ')
        # st = TC.direct_calculation_of_tc(new_gate_list, cut_list)
        st,sto = TC.direct_calculation_of_tc_queue(new_gate_list, cut_list)
        ST.append(st)
        STO.append(sto)
        print('当前最小：', minio,' ',mini ,end=' ')
        if i%100 == 0:
            minio = min(STO)
            mini = min(ST)
    print('当前最小：', min(STO), ' ', min(ST))
    return ST



'''K分区下计算最小传输代价'''
def k_count_min_st_num(gate_list,cut_list,line_sequence):
    initial_line_sequence = Initial_line_sequence(gate_list)
    new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence, gate_list, line_sequence)
    st = TC.direct_calculation_of_tc_distributed(new_gate_list, cut_list,1)
    # st = TC.direct_calculation_of_tc_queue(new_gate_list, cut_list)
    return st


def generate_partitions(remaining_letters, partition_sizes, current_partition, all_partitions):
    if len(partition_sizes) == 0:
        if len(remaining_letters) == 0:
            line = "".join(["".join(row) for row in current_partition.copy()])
            print(line,end=' ')
            st = k_count_min_st_num(gate_list,cut_list,line)
            all_partitions.append(st)
            print('当前最低传输代价：',end=' ')
            print(min(all_partitions))
        return

    current_size = partition_sizes[0]
    for comb in itertools.combinations(remaining_letters, current_size):
        next_partition = current_partition.copy()
        next_partition.append(list(comb))
        next_remaining_letters = [letter for letter in remaining_letters if letter not in comb]

        generate_partitions(next_remaining_letters, partition_sizes[1:], next_partition, all_partitions)


def generate_line2(partition_sizes):
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
            'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']

    letters = str2[0:sum(partition_sizes)]

    all_partitions = []

    generate_partitions(letters, partition_sizes, [], all_partitions)

    return all_partitions

def split_into_k_parts(N, k):
    # 计算每一份的基础大小
    base_size = N // k
    remainder = N % k  # 计算余数

    result = []
    for i in range(k):
        # 将余数分配到前 remainder 份中
        size = base_size + (1 if i < remainder else 0)
        result.append(size)

    return result


# 主函数
if __name__ == '__main__':

    # 读取qasm
    # input_filename = 'E:/python/Distributed_quantum_circuits_scheduling/reals/new_real/8bitadder.real'
    #
    # gate_list = real_to_cnot_list(input_filename)
    # print(gate_list)
    # 读取real
    # input_filename = 'E:/python/Distributed_quantum_circuits_scheduling/qasm/8bitadder2.qasm'
    # gate_list = nct_to_ncv(input_filename)

    # gate_list = [[1, 0], [2, 0], [3, 0], [4, 0], [5, 0], [6, 0], [7, 0], [8, 0], [9, 0], [10, 0], [11, 0], [12, 0], [13, 0], [14, 0], [15, 0], [16, 0], [17, 0], [18, 0], [19, 0], [20, 0], [21, 0], [22, 0], [23, 0], [24, 0], [25, 0], [26, 0], [27, 0], [28, 0], [29, 0], [30, 0], [31, 0], [2, 1], [3, 1], [4, 1], [5, 1], [6, 1], [7, 1], [8, 1], [9, 1], [10, 1], [11, 1], [12, 1], [13, 1], [14, 1], [15, 1], [16, 1], [17, 1], [18, 1], [19, 1], [20, 1], [21, 1], [22, 1], [23, 1], [24, 1], [25, 1], [26, 1], [27, 1], [28, 1], [29, 1], [30, 1], [31, 1], [3, 2], [4, 2], [5, 2], [6, 2], [7, 2], [8, 2], [9, 2], [10, 2], [11, 2], [12, 2], [13, 2], [14, 2], [15, 2], [16, 2], [17, 2], [18, 2], [19, 2], [20, 2], [21, 2], [22, 2], [23, 2], [24, 2], [25, 2], [26, 2], [27, 2], [28, 2], [29, 2], [30, 2], [31, 2], [4, 3], [5, 3], [6, 3], [7, 3], [8, 3], [9, 3], [10, 3], [11, 3], [12, 3], [13, 3], [14, 3], [15, 3], [16, 3], [17, 3], [18, 3], [19, 3], [20, 3], [21, 3], [22, 3], [23, 3], [24, 3], [25, 3], [26, 3], [27, 3], [28, 3], [29, 3], [30, 3], [31, 3], [5, 4], [6, 4], [7, 4], [8, 4], [9, 4], [10, 4], [11, 4], [12, 4], [13, 4], [14, 4], [15, 4], [16, 4], [17, 4], [18, 4], [19, 4], [20, 4], [21, 4], [22, 4], [23, 4], [24, 4], [25, 4], [26, 4], [27, 4], [28, 4], [29, 4], [30, 4], [31, 4], [6, 5], [7, 5], [8, 5], [9, 5], [10, 5], [11, 5], [12, 5], [13, 5], [14, 5], [15, 5], [16, 5], [17, 5], [18, 5], [19, 5], [20, 5], [21, 5], [22, 5], [23, 5], [24, 5], [25, 5], [26, 5], [27, 5], [28, 5], [29, 5], [30, 5], [31, 5], [7, 6], [8, 6], [9, 6], [10, 6], [11, 6], [12, 6], [13, 6], [14, 6], [15, 6], [16, 6], [17, 6], [18, 6], [19, 6], [20, 6], [21, 6], [22, 6], [23, 6], [24, 6], [25, 6], [26, 6], [27, 6], [28, 6], [29, 6], [30, 6], [31, 6], [8, 7], [9, 7], [10, 7], [11, 7], [12, 7], [13, 7], [14, 7], [15, 7], [16, 7], [17, 7], [18, 7], [19, 7], [20, 7], [21, 7], [22, 7], [23, 7], [24, 7], [25, 7], [26, 7], [27, 7], [28, 7], [29, 7], [30, 7], [31, 7], [9, 8], [10, 8], [11, 8], [12, 8], [13, 8], [14, 8], [15, 8], [16, 8], [17, 8], [18, 8], [19, 8], [20, 8], [21, 8], [22, 8], [23, 8], [24, 8], [25, 8], [26, 8], [27, 8], [28, 8], [29, 8], [30, 8], [31, 8], [10, 9], [11, 9], [12, 9], [13, 9], [14, 9], [15, 9], [16, 9], [17, 9], [18, 9], [19, 9], [20, 9], [21, 9], [22, 9], [23, 9], [24, 9], [25, 9], [26, 9], [27, 9], [28, 9], [29, 9], [30, 9], [31, 9], [11, 10], [12, 10], [13, 10], [14, 10], [15, 10], [16, 10], [17, 10], [18, 10], [19, 10], [20, 10], [21, 10], [22, 10], [23, 10], [24, 10], [25, 10], [26, 10], [27, 10], [28, 10], [29, 10], [30, 10], [31, 10], [12, 11], [13, 11], [14, 11], [15, 11], [16, 11], [17, 11], [18, 11], [19, 11], [20, 11], [21, 11], [22, 11], [23, 11], [24, 11], [25, 11], [26, 11], [27, 11], [28, 11], [29, 11], [30, 11], [31, 11], [13, 12], [14, 12], [15, 12], [16, 12], [17, 12], [18, 12], [19, 12], [20, 12], [21, 12], [22, 12], [23, 12], [24, 12], [25, 12], [26, 12], [27, 12], [28, 12], [29, 12], [30, 12], [31, 12], [14, 13], [15, 13], [16, 13], [17, 13], [18, 13], [19, 13], [20, 13], [21, 13], [22, 13], [23, 13], [24, 13], [25, 13], [26, 13], [27, 13], [28, 13], [29, 13], [30, 13], [31, 13], [15, 14], [16, 14], [17, 14], [18, 14], [19, 14], [20, 14], [21, 14], [22, 14], [23, 14], [24, 14], [25, 14], [26, 14], [27, 14], [28, 14], [29, 14], [30, 14], [31, 14], [16, 15], [17, 15], [18, 15], [19, 15], [20, 15], [21, 15], [22, 15], [23, 15], [24, 15], [25, 15], [26, 15], [27, 15], [28, 15], [29, 15], [30, 15], [31, 15], [17, 16], [18, 16], [19, 16], [20, 16], [21, 16], [22, 16], [23, 16], [24, 16], [25, 16], [26, 16], [27, 16], [28, 16], [29, 16], [30, 16], [31, 16], [18, 17], [19, 17], [20, 17], [21, 17], [22, 17], [23, 17], [24, 17], [25, 17], [26, 17], [27, 17], [28, 17], [29, 17], [30, 17], [31, 17], [19, 18], [20, 18], [21, 18], [22, 18], [23, 18], [24, 18], [25, 18], [26, 18], [27, 18], [28, 18], [29, 18], [30, 18], [31, 18], [20, 19], [21, 19], [22, 19], [23, 19], [24, 19], [25, 19], [26, 19], [27, 19], [28, 19], [29, 19], [30, 19], [31, 19], [21, 20], [22, 20], [23, 20], [24, 20], [25, 20], [26, 20], [27, 20], [28, 20], [29, 20], [30, 20], [31, 20], [22, 21], [23, 21], [24, 21], [25, 21], [26, 21], [27, 21], [28, 21], [29, 21], [30, 21], [31, 21], [23, 22], [24, 22], [25, 22], [26, 22], [27, 22], [28, 22], [29, 22], [30, 22], [31, 22], [24, 23], [25, 23], [26, 23], [27, 23], [28, 23], [29, 23], [30, 23], [31, 23], [25, 24], [26, 24], [27, 24], [28, 24], [29, 24], [30, 24], [31, 24], [26, 25], [27, 25], [28, 25], [29, 25], [30, 25], [31, 25], [27, 26], [28, 26], [29, 26], [30, 26], [31, 26], [28, 27], [29, 27], [30, 27], [31, 27], [29, 28], [30, 28], [31, 28], [30, 29], [31, 29], [31, 30], [0, 31], [31, 0], [0, 31], [1, 30], [30, 1], [1, 30], [2, 29], [29, 2], [2, 29], [3, 28], [28, 3], [3, 28], [4, 27], [27, 4], [4, 27], [5, 26], [26, 5], [5, 26], [6, 25], [25, 6], [6, 25], [7, 24], [24, 7], [7, 24], [8, 23], [23, 8], [8, 23], [9, 22], [22, 9], [9, 22], [10, 21], [21, 10], [10, 21], [11, 20], [20, 11], [11, 20], [12, 19], [19, 12], [12, 19], [13, 18], [18, 13], [13, 18], [14, 17], [17, 14], [14, 17], [15, 16], [16, 15], [15, 16]]
    # input_filename = 'E:/python/Distributed_quantum_circuits_scheduling/qasm/qft/qft64.qasm'
    # gate_list = list_str_to_int(converter_circ_from_qasm(input_filename))
    # print(mct_list)
    # gate_list = qasm_to_cnot_list(input_filename)

    gate_list = [[1,4],[0,3],[0,5],[4,0],[3,1],[2,5],[4,1],[0,2],[0,5]]
    # gate_list = read_txt1('E:/python/Distributed_quantum_circuits_scheduling/Txt/hwb50.txt')



    time_list = []

    # print(gate_list)
    # print(len(gate_list))

    for i in range(5,15):
        start_time = time.time()
        gate_list.append([0, i])
        # 量子位
        circuit_qubit = max([max(row) for row in gate_list]) + 1
        print('量子位数：' + str(circuit_qubit))

        # print('初始全局门数：' +str(len(is_global_gate(gate_list,cut_list)[1])))

        # cut_list = [12,16,12,16]
        cut_list = split_into_k_parts(circuit_qubit, 4)
        # cut_list = [12,12]

        # tc = TC.direct_calculation_of_tc(gate_list, cut_list)
        # print(tc)
        #
        # min_gg_num,min_gate_list = count_min_global_gate_num(gate_list, cut_list)
        # print('最小全局门数：' +str(min_gg_num))
        # print(min_gate_list)
        #
        # tc = TC.direct_calculation_of_tc(min_gate_list, cut_list)
        # print(tc)

        st_list = generate_line2(cut_list)

        # line_sequence_list = random_line(cut_list,10000)
        # line_sequence_list = generate_line(cut_list)
        # st = k_count_min_st_num2(gate_list,cut_list,line_sequence_list)
        # print(min(st))

        # tc = TC.direct_calculation_of_tc(gate_list, cut_list)
        # print(tc)

        # k_count_min_st_num(gate_list, cut_list, line_sequence_list)

        # print(st_list)
        # print(min(st_list))
        print(cut_list)

        end_time = time.time()
        execution_time = end_time - start_time

        print("Execution Time:", execution_time, "seconds")
        time_list.append(execution_time)
    print(time_list)
