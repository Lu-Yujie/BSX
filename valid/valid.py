import sys
import os
import glob
from subprocess import Popen, PIPE


def generate_args(binary, *params):
    arguments = [binary]
    arguments.extend(list(params))
    return arguments


def execute_binary(args):
    process = Popen(' '.join(args), shell=True, stdout=PIPE, stderr=PIPE)
    (std_output, std_error) = process.communicate()
    process.wait()
    rc = process.returncode

    # In Python 3, we need to decode the byte stream to string
    std_output = std_output.decode('utf-8')
    std_error = std_error.decode('utf-8')

    return rc, std_output, std_error


def check_correctness(binary_path, data_graph_path, query_folder_path, expected_results):
    # find all query graphs.
    query_graph_path_list = glob.glob('{0}/*.graph'.format(query_folder_path))
    query_graph_path_list = sorted(query_graph_path_list)

    # check the correctness of BSX engine.
    for query_graph_path in query_graph_path_list:
        execution_args = generate_args(binary_path, '-d', data_graph_path, '-q', query_graph_path, '-filter', 'LDF',
                                       '-order', 'GQL', '-engine', 'BSX')

        (rc, std_output, std_error) = execute_binary(execution_args)
        query_graph_name = os.path.splitext(os.path.basename(query_graph_path))[0]
        expected_embedding_num = expected_results[query_graph_name]

        if rc == 0:
            embedding_num = 0
            std_output_list = std_output.split('\n')
            for line in std_output_list:
                if '#Embeddings' in line:
                    embedding_num = int(line.split(':')[1].strip())
                    break
            if embedding_num != expected_embedding_num:
                print('BSX engine {0} is wrong. Expected {1}, Output {2}'.format(query_graph_name,
                        expected_embedding_num, embedding_num))
                # print(std_output)
                # print(std_error)
                # exit(-1)
        else:
            print('BSX engine {0} error.'.format(query_graph_name))
            print(std_output)
            print(std_error)
            exit(-1)
    print('BSX engine passed the correctness check.')


if __name__ == '__main__':
    input_binary_path = sys.argv[1]
    if not os.path.isfile(input_binary_path):
        print('The binary {0} does not exist.'.format(input_binary_path))
        exit(-1)

    # load expected results.
    input_expected_results = {}
    input_expected_results_file = 'expected_output.res'
    with open(input_expected_results_file, 'r') as f:
        for line in f:
            if line:
                result_item = line.split(':')
                input_expected_results[result_item[0].strip()] = int(result_item[1].strip())

    dir_path = os.path.dirname(os.path.realpath(__file__))
    input_data_graph_path = '{0}/data_graph/HPRD.graph'.format(dir_path)
    input_query_graph_folder_path = '{0}/query_graph/'.format(dir_path)

    check_correctness(input_binary_path, input_data_graph_path, input_query_graph_folder_path, input_expected_results)
