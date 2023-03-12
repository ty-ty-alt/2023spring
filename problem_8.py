from collections import defaultdict
from tqdm import tqdm
from count_freqs import trigram_feature_counter, trigram_feature_reader, simple_conll_corpus_iterator, sentence_iterator
from count_freqs import ViterbiTagger
import util
"""
Count n-gram frequencies in a data file and write counts to
stdout. 
"""

with open("../Code/gene.counts", "r") as H_file:
    wordlines = H_file.readlines()


def HMMemission(I_GENE, NO_GENE):
    I_G_list = list(I_GENE.keys())
    NO_G_list = list(NO_GENE.keys())

    I_only_list = [word for word in I_G_list if word not in NO_G_list]
    NO_only_list = [word for word in NO_G_list if word not in I_G_list]
    Both_list = [word for word in I_G_list if word in NO_G_list]
    assert len(Both_list) + len(I_only_list) == len(I_G_list)
    assert len(Both_list) + len(NO_only_list) == len(NO_G_list)

    I_G_total_list = [int(times) for times, tagger in I_GENE.values()]
    I_G_total = sum(I_G_total_list)

    NO_G_total_list = [int(times) for times, tagger in NO_GENE.values()]
    NO_G_total = sum(NO_G_total_list)

    HMMemissions = []

    for wordlines in I_only_list:
        HMMemissions.append({wordlines: (int(I_GENE[wordlines][0]) / I_G_total, 0)})
    for wordlines in NO_only_list:
        HMMemissions.append({wordlines: (0, int(NO_GENE[wordlines][0]) / NO_G_total)})
    for wordlines in Both_list:
        HMMemissions.append({wordlines: (int(I_GENE[wordlines][0]) / I_G_total, int(NO_GENE[wordlines][0]) / NO_G_total)})

    return HMMemissions

def rare_count_emisssion():
    # result sample: {'Ag': (0.1, 0.2)}
    with open("rare_gene.counts", "r") as rare_count_f:
        wordlines = rare_count_f.readlines()

    I_GENE = {}
    NO_GENE = {}

    for line in wordlines:
        l = line.strip()
        if l:
            parts = l.split(" ")
            if parts[1] == "WORDTAG":
                if parts[2] == "O":
                    NO_GENE.update({parts[3]: (parts[0], parts[2])})
                if parts[2] == "I-GENE":
                    I_GENE.update({parts[3]: (parts[0], parts[2])})

    result = HMMemission(I_GENE, NO_GENE)

    return result


def To_rare_dict(wordlines):
    word_dict = {}
    for l in wordlines:
        line = l.strip()
        if line:
            fields = line.split(" ")
            if fields[1] == "WORDTAG":
                count = int(fields[0])
                if count < 5:
                    word = fields[3]
                    tag = fields[2]
                    if word in word_dict:
                        current_count, current_tag = word_dict[word]
                        if count < current_count:
                            word_dict[word] = (count, tag)
                    else:
                        word_dict[word] = (count, tag)
    return word_dict



def Is_word_correct(word_dict, wordlines):
    incorrect_dict = {}
    for l in wordlines:
        line = l.strip()
        if line:
            fields = line.split(" ")
            if fields[1] == "WORDTAG":
                word = fields[3]
                count = int(fields[0])
                tag = fields[2]
                if word in word_dict:
                    stored_count, stored_tag = word_dict[word]
                    if (count, tag) != (stored_count, stored_tag):
                        total_count = count + stored_count
                        if total_count >= 5:
                            incorrect_dict[word] = total_count
    return incorrect_dict



def Toget_final_rare(write_out=False):
    word_dict = To_rare_dict(wordlines)
    incorrect_dict = Is_word_correct(word_dict, wordlines)

    final_rare_word_set = set(word_dict.keys()) - set(incorrect_dict.keys())
    final_rare_word_list = list(final_rare_word_set)

    print(len(final_rare_word_list))
    print(final_rare_word_list)

    if write_out:
        with open("rare_word_list.txt", "w") as write_file:
            for word in final_rare_word_list:
                write_file.write(str(word))
                write_file.write("###---###")
    return final_rare_word_list


def compare_tuple(a_tuple):
    if a_tuple[0] > a_tuple[1]:
        temp = "I-GENE"
    else:
        temp = "O"
    return temp


def Baseline():
    result = rare_count_emisssion()
    result_dict = {}
    for n in result:
        result_dict.update(n)

    print(result_dict)
    print(len(result_dict), type(result_dict))

    with open("../Code/gene_test.p1.out", "w") as gene_test_p1_out:

        with open("../Code/gene.test", "r") as test_file:
            wordlines = test_file.readlines()

        for line in tqdm(wordlines):
            word = line.strip()
            if word:
                if word in result_dict.keys():
                    temp = compare_tuple(result_dict[word])
                    gene_test_p1_out.write(f"{word} {temp}\n")
                else:
                    temp = compare_tuple(result_dict["_RARE_"])
                    gene_test_p1_out.write(f"{word} {temp}\n")
            else:
                gene_test_p1_out.write("\n")

Baseline()