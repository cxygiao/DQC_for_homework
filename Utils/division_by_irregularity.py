'''
不规则划分
'''

# 主函数
import itertools
import math
import random


from Utils.generate_partitions import generate_line, random_line
from Utils.min_global_gate_num import Initial_line_sequence, k_change_gate_list_by_line_sequence, \
     is_global_gate, change_gate_list_by_line_sequence
from Utils.read_qasm import converter_circ_from_qasm, list_str_to_int
from Utils.transmission_cost_calculation_more import statistics_gate_labels, transfer_qubit_list_by_gate_list, \
    count_transfer_queue, direct_calculation_of_tc




def k_count_min_gg_num(gate_list,cut_list,cut_point,line):
    initial_line_sequence = Initial_line_sequence(gate_list)
    new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence, gate_list, line)
    cuted_gate_list = []
    for i in range(len(cut_point)-1):
        cuted_gate_list.append(new_gate_list[cut_point[i]:cut_point[i+1]])
    # print(cuted_gate_list)
    gg_num = 0
    for j in range(len(cuted_gate_list)):
        gg_num_j = count_gg_num_by_line_sequence(initial_line_sequence, cut_list[j], cuted_gate_list[j], line)
        print(gg_num_j,end=' ')
        gg_num += gg_num_j
    return gg_num


'''根据线序计算全局门数'''
def count_gg_num_by_line_sequence(initial_line_sequence,cut_list, gate_list,single_line_sequence):
    new_gate_list = change_gate_list_by_line_sequence(initial_line_sequence,gate_list,single_line_sequence)
    gg_num = len(is_global_gate(new_gate_list,cut_list)[1])
    return gg_num


'''K分区下计算最小传输代价'''
def k_count_min_st_num_ahead_by_line(gate_list,cut_list,line_sequence):
    initial_line_sequence = Initial_line_sequence(gate_list)
    new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence, gate_list, line_sequence)
    st = direct_calculation_of_tc_look_ahead(new_gate_list,cut_list)
    # st = TC.direct_calculation_of_tc_queue(new_gate_list, cut_list)
    return st

#
def generate_partitions(remaining_letters, partition_sizes, current_partition, all_partitions):
    if len(partition_sizes) == 0:
        if len(remaining_letters) == 0:
            line = "".join(["".join(row) for row in current_partition.copy()])
            print(line,end=' ')
            st = k_count_min_st_num_ahead_by_line(gate_list,cut_list,line)
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


# def min_st_num(gate_list,cut_list,look_ahead_depth):
#     is_global_gate_list = is_global_gate(gate_list, cut_list)[0]  # [0,0,0,1,1,0,1,0...]
#     gg_gate_list = []
#     for i in range(len(gate_list)):
#         gg_list = []
#         if is_global_gate_list[i] == 1:
#             gg_list.append(gate_list[i])
#             for j in range(i+1, look_ahead_depth):
#                  same_qubits = list(set(gate_list[i])&set(gate_list[j]))
#                  if is_global_gate_list[j] == 1 and len(same_qubits) > 0:
#                      if len(same_qubits) == 1: #只有一个相同的量子位


#计算前瞻深度
def count_look_ahead_depth(gate_list,gate_i,statistics_gate_labels_list):
    j = gate_i+1
    while j < len(gate_list):
        if (statistics_gate_labels_list[j][0] == statistics_gate_labels_list[j][1]) and (len(list(set(gate_list[j]) & set(gate_list[gate_i]))) != 0):  # 是局部门
            break
        else:
            j+=1
    look_ahead_depth = j-gate_i-1
    # print(look_ahead_depth)
    return look_ahead_depth

#计算F
def calcluate_F(gate_list,look_ahead_depth,i,statistics_gate_labels_list,transfer_qubit):
    transfer_qubit = [transfer_qubit]
    # print(transfer_qubit)
    s_list = [] # 有没有相同量子位 有：1 ，没有：0
    e_list = [] # -1负影响：作用在传输量子位的局部门 1正影响：作用在传输量子位的全局门 0无影响：其他情况
    F = []
    for j in range(i+1,i+look_ahead_depth+1):
        if len(list(set(gate_list[j]) & set(gate_list[i]))) != 0: # 有：1
            s_list.append(1)
            if (statistics_gate_labels_list[j][0] == statistics_gate_labels_list[j][1]) and (len(list(set(gate_list[j]) & set(transfer_qubit))) != 0): # 局部门 -1
                e_list.append(-1)
            if (statistics_gate_labels_list[j][0] != statistics_gate_labels_list[j][1]) and (len(list(set(gate_list[j]) & set(transfer_qubit))) != 0): # 全局门 1
                e_list.append(1)
            else:
                e_list.append(0)
        else:
            e_list.append(0)
            s_list.append(0)

    # print(s_list)
    # print(e_list)
    for k in range(look_ahead_depth):
        F.append((look_ahead_depth-k+1)*s_list[k]*e_list[k])
    # print(F)
    return sum(F)



# 生成前瞻的transfer_qubit_list 初始量子位传输列表
def transfer_qubit_list_by_gate_list_based_look_ahead(gate_list, statistics_gate_labels_list):
    transfer_qubit_list = []
    # circuit_qubit = max([max(row) for row in gate_list]) + 1
    i = 0
    while i < len(gate_list):
        if i == len(gate_list)-1:
            transfer_qubit_list.append(gate_list[i][0])
            break
        if statistics_gate_labels_list[i][0] != statistics_gate_labels_list[i][1]:  # 是全局门
            look_ahead_depth = count_look_ahead_depth(gate_list,i,statistics_gate_labels_list)
            transfer_qubit_0 = gate_list[i][0]
            transfer_qubit_1 = gate_list[i][1]
            if calcluate_F(gate_list,look_ahead_depth,i,statistics_gate_labels_list,transfer_qubit_0) > calcluate_F(gate_list,look_ahead_depth,i,statistics_gate_labels_list,transfer_qubit_1):
                transfer_qubit_list.append(gate_list[i][0])
                i += 1
            else:
                transfer_qubit_list.append(gate_list[i][1])
                i += 1
        else: # 是局部门
            i = i + 1
    for j in range(len(transfer_qubit_list)):
        transfer_qubit_list[j] = 'q'+str(transfer_qubit_list[j])
    # print(transfer_qubit_list)
    return transfer_qubit_list


# 生成前瞻的transfer_qubit_list 初始量子位传输列表
def transfer_qubit_list_by_gate_list_based_look_ahead2(gate_list, statistics_gate_labels_list):
    transfer_qubit_list = []
    # circuit_qubit = max([max(row) for row in gate_list]) + 1
    i = 0
    while i < len(gate_list):
        if i == len(gate_list)-1:
            transfer_qubit_list.append(gate_list[i][0])
            break
        if statistics_gate_labels_list[i][0] != statistics_gate_labels_list[i][1]:  # 是全局门
            look_ahead_depth = count_look_ahead_depth(gate_list,i,statistics_gate_labels_list)
            transfer_qubit_0 = gate_list[i][0]
            transfer_qubit_1 = gate_list[i][1]
            if calcluate_F(gate_list,look_ahead_depth,i,statistics_gate_labels_list,transfer_qubit_0) >= calcluate_F(gate_list,look_ahead_depth,i,statistics_gate_labels_list,transfer_qubit_1):
                transfer_qubit_list.append(gate_list[i][0])
                i += 1
            else:
                transfer_qubit_list.append(gate_list[i][1])
                i += 1
        else: # 是局部门
            i = i + 1
    for j in range(len(transfer_qubit_list)):
        transfer_qubit_list[j] = 'q'+str(transfer_qubit_list[j])
    # print(transfer_qubit_list)
    return transfer_qubit_list


'''直接搜索'''
def direct_calculation_of_tc_look_ahead(gate_list,cut_list):

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

    return min(st0,st1)

# '''直接搜索'''
# def direct_calculation_of_tc_look_ahead(gate_list,cut_list):
#
#     # 全局门标签
#     statistics_gate_labels_list = statistics_gate_labels(gate_list, cut_list)
#     # 量子位传输列表 ['q0', 'q0', 'q2', 'q0', 'q4', 'q1', 'q0', 'q4']
#     initial_transfer_qubit_list = transfer_qubit_list_by_gate_list_based_look_ahead(gate_list, statistics_gate_labels_list)
#     # 统计合并传输队列
#     #print('######################################################################################################')
#     initial_transfer_queue = count_transfer_queue(gate_list, cut_list, initial_transfer_qubit_list,
#                                                   statistics_gate_labels_list)
#
#     print('初始传输代价：' + str(len(initial_transfer_queue) * 2),end=' ')
#     # print('######################################################################################################')
#
#     return len(initial_transfer_queue) * 2


def gate_list_to_two(gate_list,n):
    gate_list_list = []
    gate_list_list.append(gate_list[0:n])
    gate_list_list.append(gate_list[n:len(gate_list)])
    return gate_list_list

def gate_list_to_three(gate_list,n1,n2):
    gate_list_list = []
    gate_list_list.append(gate_list[0:n1])
    gate_list_list.append(gate_list[n1:n2])
    gate_list_list.append(gate_list[n2:len(gate_list)])
    return gate_list_list


'''折线'''
def direct_calculation_of_tc_look_ahead_broken(gate_list, cut_list_list):
    ST = []
    for l in range(len(gate_list)-1):
        gate_list_list = gate_list_to_two(gate_list,l+1)
        # print(gate_list_list)
        St = []
        for i in range(len(gate_list_list)):
            try:
                # 全局门标签
                statistics_gate_labels_list = statistics_gate_labels(gate_list_list[i], cut_list_list[i])
                # 量子位传输列表 ['q0', 'q0', 'q2', 'q0', 'q4', 'q1', 'q0', 'q4']
                initial_transfer_qubit_list = transfer_qubit_list_by_gate_list_based_look_ahead(gate_list_list[i],
                                                                                                statistics_gate_labels_list)
                # 统计合并传输队列
                # print('######################################################################################################')
                initial_transfer_queue = count_transfer_queue(gate_list_list[i], cut_list_list[i],
                                                              initial_transfer_qubit_list,
                                                              statistics_gate_labels_list)

                st = len(initial_transfer_queue) * 2
                St.append(st)
            except Exception as e:
                pass

        ST.append(sum(St))
    print(ST)
    return min(ST)


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


'''K分区下计算最小传输代价 前瞻'''
def k_count_min_st_num_ahead(gate_list,cut_list,line_sequence_list):
    initial_line_sequence = Initial_line_sequence(gate_list)
    ST = []
    lenl = len(line_sequence_list)
    print(lenl)
    mini = len(gate_list)*2
    for i in range(lenl):
        print(i, end=':')
        print(line_sequence_list[i],end='')
        new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,line_sequence_list[i])
        st = direct_calculation_of_tc_look_ahead(new_gate_list, cut_list)
        ST.append(st)
        print('当前最小：',mini )
        if i%100 == 0:
            mini = min(ST)

    print('最小：', min(ST))
    min_st_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,line_sequence_list[ST.index(min(ST))])
    # 绘制电路
    # draw_color_circuit(min_st_gate_list)
    return ST




'''K分区下计算最小传输代价 前瞻'''
def k_count_min_st_num_ahead_by_broken(gate_list,cut_list_list,line_sequence_list):
    initial_line_sequence = Initial_line_sequence(gate_list)
    ST = []
    lenl = len(line_sequence_list)
    print(lenl)
    mini = len(gate_list)*2
    for i in range(lenl):
        print(i, end=':')
        print(line_sequence_list[i],end='')
        new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,line_sequence_list[i])
        st = direct_calculation_of_tc_look_ahead_broken(new_gate_list, cut_list_list)
        ST.append(st)
        print('当前最小：',mini )
        if i%100 == 0:
            mini = min(ST)

    print('最小：', min(ST))
    min_st_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,line_sequence_list[ST.index(min(ST))])
    # 绘制电路
    # draw_color_circuit(min_st_gate_list)
    return ST




# 生成染色体
def generate_chromosome(length, non_zero_genes):
    # 负载均衡系数
    k = math.ceil(qubit/Part_num*Load_balancing)
    print(k)
    # 生成长度为length的染色体
    chromosome = [0] * length
    while 1:
        # 随机选择non_zero_genes个基因位置，并赋予随机值
        non_zero_indices = random.sample(range(length), non_zero_genes)
        for index in non_zero_indices:
            chromosome[index] = random.randint(-k, k)  # 假设基因的范围为-10到10

        if sum(chromosome) <= k and sum(chromosome) >= -k:
            return chromosome



#生成个体
def generate_individual(length, non_zero_genes):
    individual = []
    for i in range(Part_num-1):
        individual.append(generate_chromosome(length, non_zero_genes))
    for j in range(len(individual)):
        cumulative_list = [sum(individual[j][:i+1])+cut_list[j] for i in range(len(individual[j]))]

    return individual


#生成种群
def generate_community(length, non_zero_genes):
    community = []
    for i in range(Community_size):
        community.append(generate_individual(length, non_zero_genes))
    return community

# 个体适应度函数 st+gg
def fitness_function(gate_list,individual_lsit,cut_list):
    st = 0
    partation_list = []
    for individual in individual_lsit:
        st += sum(abs(num) for num in individual if num != 0)
    # for i in range(len(gate_list)):
    #      if sum(abs(individual_lsit[j][i]) for j in range(len(individual_lsit))) == 0:
    #          partation_list.append(cut_list)
    #      else:
    #          partation_list.append(cut_list)
    print(individual_lsit)
    for j in range(len(individual_lsit)):
        cumulative_list = [sum(individual_lsit[j][:i+1])+cut_list[j] for i in range(len(individual_lsit[j]))]
        print(cumulative_list)
    return st

# 三分区随机切割
def split_number_3(n, p):
    # 将 n 分成三个数
    a = random.randint(0, n)
    b = random.randint(0, n - a)
    c = n - a - b

    # 调整三个数，使得它们的差异最大不超过 p
    while max(a, b, c) - min(a, b, c) != p:
        a = random.randint(0, n)
        b = random.randint(0, n - a)
        c = n - a - b

    return [a,b,c]

# 四分区随机切割
def split_number_4(n, p):
    # 将 n 分成四个数
    a = random.randint(0, n)
    b = random.randint(0, n - a)
    c = random.randint(0, n - a - b)
    d = n - a - b - c

    # 调整四个数，使得它们的差异最大不超过 p
    while max(a, b, c, d) - min(a, b, c, d) != p:
        a = random.randint(0, n)
        b = random.randint(0, n - a)
        c = random.randint(0, n - a - b)
        d = n - a - b - c

    return [a, b, c, d]

# K分区随机切割
def split_number_k(n, k, p):
    # 将 n 分成 k 个数
    parts = [random.randint(0, n) for _ in range(k - 1)]
    parts.append(n - sum(parts))

    # 调整 k 个数，使得它们的差异最大不超过 p
    while max(parts) - min(parts) != p:
        parts = [random.randint(0, n) for _ in range(k - 1)]
        parts.append(n - sum(parts))

    return parts

def new_cut_list_list(n, p, random_num):
    cut_list_list = []
    for i in range(random_num):
        cut_list_list.append(split_number_4(n, p))
    return cut_list_list


def random_line_ahead(random_num):
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
            'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    letters = str2[0:qubit]
    all_partitions = []
    for i in range(random_num):
        random.shuffle(letters)
        all_partitions.append(''.join(letters))
    return all_partitions


'''K分区下计算最小传输代价 前瞻'''
def k_count_min_st_num_ahead_by_line_sequence_list(gate_list,cut_list_list,random_num):
    initial_line_sequence = Initial_line_sequence(gate_list)
    ST = []
    line_sequence_list = random_line_ahead(random_num)
    lenl = len(line_sequence_list)
    print(lenl)
    mini = len(gate_list)*2
    for i in range(lenl):
        print(i, end=':')
        print(line_sequence_list[i],end='')
        print(cut_list_list[i],end='')
        new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,line_sequence_list[i])
        st = direct_calculation_of_tc_look_ahead(new_gate_list, cut_list_list[i])
        ST.append(st)
        print('当前最小：',mini )
        if i%100 == 0:
            mini = min(ST)

    print('最小：', min(ST))
    min_st_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,line_sequence_list[ST.index(min(ST))])
    # 绘制电路
    # draw_color_circuit(min_st_gate_list)
    return ST



if __name__ == '__main__':

    # 参数定义
    Load_balancing = 0.1
    Part_num = 3
    Community_size = 10
    st_weight = 0.7
    gg_weight = 0.3
    Number_of_partitions = 3

    # input_filename = 'E:/python/Distributed_quantum_circuits_scheduling/reals/rd53_138.real'
    # gate_list = real_to_cnot_list(input_filename)
    # gate_list = read_txt1('E:/python/Distributed_quantum_circuits_scheduling/Txt/hwb50.txt')

    input_filename = 'E:/python/Distributed_quantum_circuits_scheduling/qasm/CZ_circuit/ham7_299_fj_opt_czopt.qasm'
    gate_list =  list_str_to_int(converter_circ_from_qasm(input_filename))

    qubit = max([max(row) for row in gate_list]) + 1
    print(gate_list)
    print('量子门数：',len(gate_list))
    print('量子位数：'+str(qubit))


    cut_list = split_into_k_parts(qubit, 2)
    # cut_list = split_number(qubit, 3)
    print(cut_list)



    # individual_lsit = generate_individual(len(gate_list), 3)
    #
    # fitness_function(gate_list, individual_lsit, cut_list)

    # 遍历排列组合
    line_sequence_list = generate_line(cut_list)

    # cut_list_list = new_cut_list_list(qubit, 2, 1000000)
    # k_count_min_st_num_ahead_by_line_sequence_list(gate_list,cut_list_list,1000000)

    # line_sequence_list = random_line_ahead(cut_list, 100000)
    #
    st = k_count_min_st_num_ahead(gate_list,cut_list,line_sequence_list)


    # cut_list_list1 = [[3, 2], [2, 3]]
    # cut_list_list2 = [[2, 3], [3, 2]]
    # # # direct_calculation_of_tc_look_ahead_broken(gate_list, [[2, 3], [3, 2]])
    # #
    # # # 遍历排列组合
    # line_sequence_list = generate_line(cut_list)
    # st1 = k_count_min_st_num_ahead_by_broken(gate_list,cut_list_list1,line_sequence_list)
    # st2 = k_count_min_st_num_ahead_by_broken(gate_list, cut_list_list2, line_sequence_list)
    # print(min(st1)+ 1)
    # print(min(st2) + 1)

    # st_list = generate_line2(cut_list)

    # cut_list = [[7,5,5],[5,5,7],[5,6,6]]
    # cut_point = [0,5,9,len(gate_list)]

    # if len(cut_list) == len(cut_point)-1:
    #     generate_line2(cut_list[0])


    print(input_filename)






