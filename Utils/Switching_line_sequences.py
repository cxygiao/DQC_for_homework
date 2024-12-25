'''
功能：迭代每种线序下的传输代价，输出最小传输代价
开始时间：2024.12.05
作者：陈新宇
'''

import datetime
import re
from multiprocessing.pool import Pool
import itertools
import random
import Utils.transmission_cost_calculation_more as TC
from Utils.read_qasm import converter_circ_from_qasm, count_num_of_qubit


'''
将gate_list全部转换为int
'''
def list_str_to_int(gate_list):
    new_gate_list = []
    for i in range(len(gate_list)):
        son_new_gate_list = list(map(int, gate_list[i]))
        new_gate_list.append(son_new_gate_list)
    return new_gate_list



'''根据线序改变gate_list'''
def change_gate_list_by_line_sequence(initial_line_sequence,gate_list,all_line_sequence):
    new_gate_list = list_str_to_int(gate_list)
    circuit_qbit =  max([max(row) for row in new_gate_list]) + 1
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

'''
把字母变成对于的数字
'''
def letter_to_number(letter):
    number = ord(letter) - 65
    return number

'''生成初始线序'''
def new_initial_line_sequence(circuit_qubit):
    initial_line_sequence = ''
    # 最多支持 26量子位
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z']
    for i in range(int(circuit_qubit)):
        # str1：ABCDEFGHIJKL.... 共qbit个
        initial_line_sequence = initial_line_sequence + str2[i]
    return initial_line_sequence


'''随机扰动线序'''
def randomly_perturbed_line_sequences(line_sequence):
    line_sequence = list(line_sequence)
    random.shuffle(line_sequence)
    return line_sequence


# 多进程的子函数
def multiprocessing_count_TC(initial_line_sequence,gate_list,cut_list):
    single_line_sequence = randomly_perturbed_line_sequences(initial_line_sequence)
    #print(single_line_sequence, end=' ')
    gate_list = change_gate_list_by_line_sequence(initial_line_sequence, gate_list, single_line_sequence)
    # print(gate_list, end=' ')
    tc = TC.taboo_search2(gate_list, cut_list)
    # print(tc)
    return tc,gate_list


# 多进程
def multi_processing_count_TC(TC_list):
    num_process = 10
    results = []
    pool = Pool(processes=num_process)
    for n in range(num_process):
        result = pool.apply_async(multiprocessing_count_TC, (initial_line_sequence, gate_list, cut_list))
        results.append(result)
    pool.close()
    pool.join()
    res = [i.get() for i in results]
    for i in range(len(res)):
        TC_list.append(res[i])
    return TC_list


'''TC_list中找最低传输代价的gate_list'''
def search_min_tc_in_gate_list(TC_list):
    cost_list = []
    for i in range(len(TC_list)):
        cost_list.append(TC_list[i][0])
    return TC_list[cost_list.index(min(cost_list))][1],min(cost_list)



def look_ahead_by_line_switching(gate_list,cut_list):
    from Utils.look_ahead_and_schedule import direct_calculation_of_tc_look_ahead_new
    # 在这里使用导入的函数
    # 初始线序
    initial_line_sequence = new_initial_line_sequence(circuit_qubit)
    # print(initial_line_sequence)

    # for i in range(3):
    #     single_line_sequence = randomly_perturbed_line_sequences(initial_line_sequence)
    #     print(single_line_sequence,end=' ')
    #     gate_list = change_gate_list_by_line_sequence(initial_line_sequence,gate_list,single_line_sequence)
    #     print(gate_list,end=' ')
    #     print('传输代价：' + str(TC.taboo_search2(gate_list, cut_list)))
    #
    # TC_list = []
    # min_transmission_cost = TC.taboo_search2(gate_list,cut_list)
    # for i in range(50):
    #     TC_list = multi_processing_count_TC(TC_list)
    #     one_min_gate_list = search_min_tc_in_gate_list(TC_list)[0]
    #     one_min_transmission_cost = search_min_tc_in_gate_list(TC_list)[1]
    #     if one_min_transmission_cost < min_transmission_cost:
    #         min_transmission_cost = one_min_transmission_cost
    #         min_gate_list = one_min_gate_list
    # print(TC_list)
    # print(str(min_gate_list)+' 传输代价：'+str(min_transmission_cost))

    min_global_gate_num = len(TC.is_global_gate(gate_list, cut_list)[1])
    min_tc = len(TC.is_global_gate(gate_list, cut_list)[1]) * 2
    items = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
             'V', 'W', 'X', 'Y', 'Z']

    # 全排列
    for p in itertools.permutations(items[0:circuit_qubit]):
        line_sequence = (''.join(p))
        print(line_sequence, end=' ')
        new_gate_list = change_gate_list_by_line_sequence(initial_line_sequence, gate_list, line_sequence)  # 门
        # 未前瞻
        # tc = TC.direct_calculation_of_tc(new_gate_list,cut_list)
        # 前瞻
        tc = direct_calculation_of_tc_look_ahead_new(new_gate_list, cut_list)
        global_gate_num = len(TC.is_global_gate(new_gate_list, cut_list)[1])
        print(' 全局门数:', global_gate_num, end=' ')
        print(' 当前传输代价:', tc, end=' ')
        # 找最低传输代价
        if tc < min_tc:
            best_gate_list = new_gate_list
            min_tc = tc
        print(' 当前最小传输代价:', min_tc)
    print(best_gate_list)
    print(min_tc)
    #     # 找最少全局门
    # if global_gate_num < min_global_gate_num:
    #         min_global_gate_num = global_gate_num
    #         best_gate_list = new_gate_list
    # print(best_gate_list)
    # print(min_global_gate_num)
    return min_tc


# 主函数
if __name__ == '__main__':
    start_time = datetime.datetime.now()

    input_filename = 'E:/python/Distributed_quantum_circuits_scheduling/qasm/mod5adder_127.qasm'
    gate_list = list_str_to_int(converter_circ_from_qasm(input_filename))
    print(gate_list)

    # print('量子位数：'+str(count_num_of_qubit(gate_list)))
    # 门列表

    # gate_list = [[0, 1], [2, 1], [1, 2], [3, 2], [2, 3], [4, 3], [3, 4], [4, 3], [4, 2], [4, 1], [4, 0], [0, 4], [0, 4], [1, 4],
    #  [1, 4], [2, 4], [2, 4], [1, 4], [1, 2]]


    # 画电路图
    # draw_color_circuit(gate_list)

    # gate_list = [[0, 1], [0, 1], [0, 2], [0, 2], [0, 3], [0, 3], [0, 6], [0, 6], [0, 7], [0, 7], [0, 8], [0, 8], [0, 9], [0, 9], [0, 10], [0, 10], [0, 11], [0, 11], [1, 2], [1, 2], [1, 3], [1, 3], [1, 6], [1, 6], [1, 7], [1, 7], [1, 8], [1, 8], [1, 9], [1, 9], [1, 10], [1, 10], [1, 11], [1, 11], [2, 3], [2, 3], [2, 6], [2, 6], [2, 7], [2, 7], [2, 8], [2, 8], [2, 9], [2, 9], [2, 10], [2, 10], [2, 11], [2, 11], [3, 6], [3, 6], [3, 7], [3, 7], [3, 8], [3, 8], [3, 9], [3, 9], [3, 10], [3, 10], [3, 11], [3, 11], [6, 7], [6, 7], [6, 8], [6, 8], [6, 9], [6, 9], [6, 10], [6, 10], [6, 11], [6, 11], [7, 8], [7, 8], [7, 9], [7, 9], [7, 10], [7, 10], [7, 11], [8, 9], [8, 9], [8, 10], [8, 10], [8, 11], [8, 11], [9, 10], [9, 10], [9, 11], [9, 11], [10, 11], [10, 11]]
    # 量子位
    circuit_qubit = max([max(row) for row in gate_list]) + 1
    # 切割点列表， [2,2,3]表示2，2，3个量子位分别为一个分区
    cut_list = [3,3]

    look_ahead_by_line_switching(gate_list, cut_list)


    end_time = datetime.datetime.now()
    print(end_time-start_time)





