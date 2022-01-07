from string import punctuation
from nltk.cluster.util import cosine_distance
from nltk.tag.perceptron import PerceptronTagger
from nltk.corpus import wordnet as wn
from nltk.tag import SequentialBackoffTagger,_pos_tag
from nltk.stem.porter import PorterStemmer
from nltk.probability import FreqDist
from nltk.corpus import stopwords
from nltk import edit_distance
from nltk import download as nltk_download
import pickle
from statistics import median,mean
from math import log10
import os
import re
from sys import platform
from pathlib import Path
import requests
import urllib.request as request
import shutil
from contextlib import closing
from gzip import open as gzip_open


if platform.startswith('win'):    SPLITTER = '\\'
else:                                SPLITTER = '/'

unifunc_folder = os.path.abspath(os.path.dirname(__file__)).split(SPLITTER)[0:-1]
unifunc_folder = SPLITTER.join(unifunc_folder)+SPLITTER


def run_unifunc(str1,str2,verbose=False,console_output=False,threshold=None):
    if threshold: threshold=float(threshold)
    nlp = UniFunc()
    nlp.get_similarity_score(str1, str2, verbose=verbose,threshold=threshold,console_output=console_output)
    if verbose:
        print_citation_unifunc()

def print_citation_unifunc():
    paper_doi='https://doi.org/10.1515/hsz-2021-0125'
    separator='#####################################################################################################################'
    res=f'{separator}\n# Thank you for using UniFunc, please make sure you cite the respective paper {paper_doi} #\n{separator}'
    print(res)

def print_version(user,project):
    import requests
    response = requests.get(f"https://api.github.com/repos/{user}/{project}/releases/latest")
    json=response.json()
    if 'name' in json:
        print(f'{project}\'s latest release is:',json['name'])
    else:
        print('No release available')

def run_example():
    nlp = UniFunc()
    print('####################################')
    print('###############TEST 1###############')
    print('####################################')
    str1 = 'Responsible for trypanothione reduction'
    str2='Protein associated with trypanothione reductase activity'
    print('Similarity score:',nlp.get_similarity_score(str1,str2,verbose=True))
    print('####################################')
    print('###############TEST ###############')
    print('####################################')
    str1='Leghemoglobin reductase activity K0002 (EC 0.0.0.0) ID12345 PRK10411.1  '
    str2='Protein associated with trypanothione reductase activity (K0001) ID6789'
    print('Similarity score:',nlp.get_similarity_score(str1,str2,verbose=True))
    print_citation_unifunc()







class Word_Weighter():

    def process_uniprot(self,line):
        if 'Entry' not in line:
            line = line.split('\t')
            header,annotation,function_description=line
            function_description = function_description.replace('FUNCTION: ', '')
            function_rule = re.search('\{.*\}\.', function_description)
            if function_rule:
                function_description = function_description.replace(function_rule.group(), '')
            return header,annotation,function_description

    def is_float(self,x):
        try:
            float(x)
            return True
        except:
            return False

    def is_abbreviation(self, x):
        if re.search('#[a-z]?[A-Z]+', x):
            return True
        return False

    def generate_n_grams(self,list_words):
        for n_range in self.n_grams_range:
            grams = [list_words[i:i + n_range] for i in range(len(list_words) - n_range + 1)]
            for g in grams:
                yield ' '.join(g)

    def is_good_word(self,word,stop_words):
        if word in stop_words or \
            not word or\
            len(word) == 1 or \
            re.match('\.?\d+', word):
            return False
        return True


    def build_frequency_dict(self):
        self.load_word_counter_pickle()
        if not self.word_counter:
            print('Building frequency counter into ', self.word_counter_pickled, flush=True)
            self.build_freq_dict_all()
            self.get_metrics_idfs()
            self.save_word_counter_pickle()


    def get_metrics_idfs(self):
        all_idfs=[]
        for word in self.word_counter:
            word_idf=self.calculate_idf(word)
            all_idfs.append(word_idf)

    def save_word_counter_pickle(self):
        with open(self.word_counter_pickled, 'wb') as handle:
            pickle.dump((self.word_counter,self.document_counter), handle)

    def load_word_counter_pickle(self):
        if os.path.exists(self.word_counter_pickled):
            with open(self.word_counter_pickled, 'rb') as handle:
                self.word_counter,self.document_counter= pickle.load(handle)

    def save_verbs_set_pickle(self):
        with open(self.verbs_set_pickled, 'wb') as handle:
            pickle.dump(self.verbs_set, handle)

    def load_verbs_set_pickle(self):
        if os.path.exists(self.verbs_set_pickled):
            with open(self.verbs_set_pickled, 'rb') as handle:
                self.verbs_set= pickle.load(handle)

    def get_verbs_set(self):
        self.load_verbs_set_pickle()
        if not self.verbs_set:
            self.verbs_set=set()
            with open(self.verbs_txt) as file:
                line=file.readline()
                while line:
                    line=line.strip('\n')
                    line=line.split('\t')
                    line=[i.strip() for i in line if i]
                    self.verbs_set.update(line)
                    line=file.readline()
            self.save_verbs_set_pickle()

    def min_max_scale_nlp(self,X,minX,maxX):
        if minX==maxX: return 1
        return (X-minX)/(maxX-minX)


    #is token too common in the corpus?
    def calculate_idf(self, word):
        lower_word = word.lower()
        if lower_word in self.words_to_remove:
            res = 0
        elif self.is_float(lower_word):
            res = 0
        elif lower_word in set(punctuation):
            res = 0
        elif lower_word in self.word_counter:
            N_docs_with_token = self.word_counter[lower_word]
            res = self.document_counter / N_docs_with_token
        elif '-' in lower_word:
            N_docs_with_token = self.default_word_count
            res = self.document_counter / N_docs_with_token
        else:
            res = 1
        return res

    #when coming from the main nlp, with all of its pre processing
    #if the word appears multiple times it should be more important?
    def calculate_tf_idf(self,list_words):
        all_words=[]
        for word in self.generate_n_grams(list_words):
            all_words.append(word)
        sentence_word_counter={}
        for word in all_words:
            if word not in sentence_word_counter: sentence_word_counter[word]=0
            sentence_word_counter[word]+=1
        res={}
        for word in sentence_word_counter:
            n_token_in_doc=sentence_word_counter[word]
            total_tokens_in_doc=len(all_words)
            tf=n_token_in_doc/total_tokens_in_doc
            idf = self.calculate_idf(word)
            res[word]=tf*idf+0.01
        return res

    #we scale the weights to understand which tokens are more important within each sentence
    def calculate_scaled_tf_idf(self,tf_idfs,min_max=False):
        #log10 because some idfs are very divergent, we want to make their distance a bit closer
        tf_idfs= {i:log10(tf_idfs[i]) for i in tf_idfs if tf_idfs[i]}
        if min_max:
            tf_idfs_vals=tf_idfs.values()
            maxX = max(tf_idfs_vals)
            minX = min(tf_idfs_vals)
            for word in tf_idfs:
                tf_idfs[word] = self.min_max_scale_nlp(tf_idfs[word], minX, maxX)
                if not tf_idfs[word]: tf_idfs[word]+=+0.01
        return tf_idfs




####PREPROCESSING

class Pre_Processer():
    def remove_ecs(self, string_to_search, required_level=3):
        matches = []
        removed_ecs = set()
        # greedy match of confounders
        ec_pattern = re.compile('\d(\.(-|\d{1,3}|([a-zA-Z]\d{1,3}))){2,3}')
        search = re.finditer(ec_pattern, string_to_search)
        for i in search:
            ec = i.group()
            passed = False
            start = i.span()[0]
            end = i.span()[1]
            if len(string_to_search) > end + 1:
                if string_to_search[start - 1] != '.' and string_to_search[end - 1] != '.' and not re.match(
                        '\.|[a-zA-Z]|\d{1,3}', string_to_search[end + 1]) and not re.match('-', string_to_search[end]):
                    passed = True
            else:
                if string_to_search[start - 1] != '.':
                    passed = True
            if passed:
                if ec.count('.') >= required_level - 1:
                    if ec.count('.') + 1 - ec.count('-') >= required_level:
                        matches.append([start, end])
                        removed_ecs.add(ec)
        removed_length = 0
        current_str = str(string_to_search)
        for m in matches:
            start, end = m
            start -= removed_length
            end -= removed_length
            current_str = current_str[:start] + current_str[end:]
            removed_length += end - start
        to_remove_pattern = re.compile(
            '(\(EC\s? \)|\(EC:\s?\)|\(ec\s?\)|\(ec:\s?\)|\[EC\s?\]|\[EC:\s?\]|\[ec\s?\]|\[ec:\s?\])')
        search = re.search(to_remove_pattern, current_str)
        if search:
            current_str = current_str.replace(search.group(), '')
        return current_str, removed_ecs

    def remove_pattern(self, string_to_search, pattern):
        #this sort of patterns should always be preceded by a not a letter or a digit. or preceded by space ( or [
        negative_pattern='[A-Za-z\-\d+\_\.\']'
        compiled_pattern=re.compile('(?<!('+negative_pattern+'))'+pattern+'(?!('+negative_pattern+'))')
        patterns_removed = set()
        search = re.search(compiled_pattern, string_to_search)
        while search:
            patterns_removed.add(search.group())
            start = search.span()[0]
            end = search.span()[1]
            #if string_to_search[start + 1] == '(': start += 2
            #if string_to_search[end - 1] == ')': end -= 1
            string_to_search = list(string_to_search)
            string_to_search[start:end] = ''
            string_to_search = ''.join(string_to_search)
            search = re.search(compiled_pattern, string_to_search)
        return string_to_search, patterns_removed

    def convert_to_arabic_digits(self, roman_digit):
        roman_numerals = [
            ('M', 1000),
            ('CM', 900),
            ('D', 500),
            ('CD', 400),
            ('C', 100),
            ('XC', 90),
            ('L', 50),
            ('XL', 40),
            ('X', 10),
            ('IX', 9),
            ('V', 5),
            ('IV', 4),
            ('I', 1)
        ]
        ix = 0
        result = 0
        while ix < len(roman_digit):
            for k, v in roman_numerals:
                if roman_digit.startswith(k, ix):
                    result += v
                    ix += len(k)
                    break
            else:
                raise ValueError('Invalid Roman number.')
        return result

    def replace_roman_numerals(self, string_to_search):
        # we wont use high roman numerals since they dont usually come up in this scenario
        roman_pattern = re.compile('[^a-zA-Z0-9][IV]+[^a-zA-Z0-9\.]')
        search = re.search(roman_pattern, string_to_search)
        while search:
            start = search.span()[0]
            end = search.span()[1]
            roman_digit = re.search('[IV]+', search.group()).group()
            string_to_search = list(string_to_search)
            converted_number = search.group().replace(roman_digit.upper(),
                                                      str(self.convert_to_arabic_digits(roman_digit)))
            string_to_search[start:end] = converted_number
            string_to_search = ''.join(string_to_search)
            search = re.search(roman_pattern, string_to_search)

        return string_to_search

    def replace_punctuation(self, sentence):
        temp_sentence = str(sentence)
        temp_sentence = temp_sentence.replace('-->', 'to')
        punctuation_set = set(punctuation)
        for i in ['\'', '-','_', '.', ',', '+', '(', ')', '[', ']','/']:
            punctuation_set.remove(i)
        punctuation_set.add(', ')
        punctuation_set.add(' - ')
        for p in punctuation_set:
            temp_sentence = temp_sentence.replace(p, ' ')
        # some terms are wrongly separated ex:GDP-4- dehydro
        temp_sentence = temp_sentence.replace('- ', '-')
        return temp_sentence

    def remove_bad_pattern_ions(self, string_to_search):
        bad_ion_pattern = re.compile('\(\d\s\)')
        search = re.search(bad_ion_pattern, string_to_search)
        while search:
            ion_str = search.group()
            new_ion_str = ion_str.replace(' ', '')
            string_to_search = string_to_search.replace(ion_str, new_ion_str)
            search = re.search(bad_ion_pattern, string_to_search)
        return string_to_search

    def simplify_ions(self, string_to_search):
        ion_pattern = re.compile('(\(\d\([\+\-]\)\))|(\(\d\)\([\+\-])\)')
        search = re.search(ion_pattern, string_to_search)
        while search:
            ion_str = search.group()
            new_ion_str = '(' + ion_str.replace('(', '').replace(')', '') + ')'
            string_to_search = string_to_search.replace(ion_str, new_ion_str)
            search = re.search(ion_pattern, string_to_search)
        return string_to_search

    def get_parentheses_above_range(self, starting_range, p2_ranges):
        res = []
        for p2 in p2_ranges:
            if p2[0] > starting_range[0]: res.append(p2)
        return res

    def get_parentheses_pairs(self, p1_ranges, p2_ranges):
        res = []
        remaining = []

        if not p1_ranges or not p2_ranges: return res, remaining
        while p1_ranges:
            current_p1 = p1_ranges.pop(-1)
            p2_above_range = self.get_parentheses_above_range(current_p1, p2_ranges)
            if not p2_above_range:
                remaining = list(p1_ranges)
                break
            p2_for_p1 = p2_above_range[0]
            p2_ranges.remove(p2_for_p1)
            res.append([current_p1, p2_for_p1])
        return res, remaining

    def remove_parentheses(self, string_to_search):
        p1_pattern = re.compile('\(')
        p2_pattern = re.compile('\)')
        p1_search = re.finditer(p1_pattern, string_to_search)
        p2_search = re.finditer(p2_pattern, string_to_search)
        p1_search = [i.span() for i in p1_search]
        p2_search = [i.span() for i in p2_search]
        res = list(string_to_search)
        # when parentheses are loose we just removed them
        if not p2_search:
            for p1 in p1_search:
                res[p1[0]] = ''
        if not p1_search:
            for p2 in p2_search:
                res[p2[0]] = ''
        p_pairs, remaining = self.get_parentheses_pairs(p1_search, p2_search)
        for r in remaining:
            res[r[0]] = ''
        to_remove = []
        for i in p_pairs:
            if string_to_search[i[0][0] - 1:i[0][1]] == ' (' and (string_to_search[i[1][0]:i[1][1] + 1] == ') '
                                                                  or string_to_search[i[1][0]:i[1][1] + 2] == '). '
                                                                  or string_to_search[i[1][0]:i[1][1] + 2] == ').\t'
                                                                  or string_to_search[i[1][0]:i[1][1] + 1] == ')\t'):
                to_remove.append(i)
        for t in to_remove:
            res[t[0][0]] = '#SPLIT#'
            res[t[1][0]] = '#SPLIT#'
        res = ''.join(res)
        return res

    def unite_terms(self, string_to_process):
        #this is to unite certain tokens e.g. flavoprotein,A to flavoprotein A
        string_list = string_to_process.split()
        res = []
        c = 0
        for i in range(len(string_list)):
            if i == 0:
                res.append(string_list[i])
            else:
                if (len(string_list[i]) == 1 and string_list[i].isupper()) or re.search('\d+\s', string_list[i]):
                    c += 1
                    res[i - c] += ' ' + string_list[i]
                else:
                    res.append(string_list[i])
        c = 0
        res=[word for word in res if word]
        for i in range(len(res)):
            if res[i]:
                if res[i][0] == '-':
                    res[i] = res[i].lstrip('-')
                    c += 1
            if res[i]:
                if res[i][-1] == '-':
                    res[i] = res[i].rstrip('-')
                    c += 1
            if res[i]:
                if (res[i].count('(') == 1 and res[i].count(')') == 0) or \
                        (res[i].count(')') == 1 and res[i].count('(') == 0):
                    c += 1
                    res[i] = res[i].strip(')')
        return res

    def get_token_to_merge(self, list_of_tokens, to_add_pos):
        for i in range(len(list_of_tokens)):
            if i > to_add_pos:
                passed = True
                for t in list_of_tokens[::-1][i]:
                    if '#' in t: passed = False
                if passed:
                    list_of_tokens[::-1][i].extend(list_of_tokens[::-1][to_add_pos][1:])
                    to_remove = list_of_tokens[::-1].pop(to_add_pos)
                    list_of_tokens.remove(to_remove)
                    return list_of_tokens

    def connect_gapped_token(self, list_of_tokens):
        to_append = []
        parentheses_tokens=[]
        for i in range(len(list_of_tokens)):
            if list_of_tokens[i][0] == '#':
                to_append.append(i)
        for i in to_append[::-1]:
            self.get_token_to_merge(list_of_tokens, len(list_of_tokens) - i - 1)
        for li in range(len(list_of_tokens)):
            is_parentheses=False
            for inner_str in list_of_tokens[li]:
                if inner_str[0]=='#': is_parentheses=True
            if len(list_of_tokens[li])==1  and self.is_abbreviation(list_of_tokens[li][0]):
                is_parentheses=False
                list_of_tokens[li][0]=list_of_tokens[li][0].strip('#')
                list_of_tokens[li-1].append(list_of_tokens[li][0])
                list_of_tokens[li]=[]
            if is_parentheses: parentheses_tokens.append(list_of_tokens[li])

        res = []
        parentheses_res = []
        for lt in list_of_tokens:
            temp = []
            if lt not in parentheses_tokens:
                for token in lt:
                    if token:
                        temp.append(token)
                if temp:
                    res.append(temp)
        for i in range(len(parentheses_tokens)):
            for j in range(len(parentheses_tokens[i])):
                parentheses_tokens[i][j]=parentheses_tokens[i][j].strip('#')
        for lt in parentheses_tokens:
            temp = []
            for token in lt:
                if token:
                    token=token.replace('\'\'','')
                    temp.append(token)
            if temp:
                parentheses_res.append(temp)
        return res,parentheses_res


    def split_by_alternatives(self,string_to_process):
        if '/' not in string_to_process: return string_to_process
        temp=str(string_to_process)
        res=[]
        alternatives=string_to_process.split('/')
        suffix=' '.join(alternatives[-1].split()[1:]).strip()
        temp=temp.replace(suffix,'').strip()
        add_split = False
        for i in temp.split('/'):
            if i + ' ' + suffix not in res:
                res.append(i + ' ' + suffix)
        for i in range(len(res)):
            if add_split:res[i]='#SPLIT#'+res[i]
            if res[i].endswith('#SPLIT#'):add_split=True
            else:add_split=False
        return '!NEWLINE!'.join(res)

    def remove_extra_parentheses(self,current_str):
        res=str(current_str)
        if current_str.count('(') != current_str.count(')'):
            res = res.strip('()')
        if current_str.count('[') != current_str.count(']'):
            res = res.strip('[]')
        res = res.strip('.')
        return res

    def split_by_parentheses(self,string_to_process):
        #this unites pharses split by tokens in parentheses
        n_parentheses=string_to_process.count('#SPLIT#')/2
        split_parentheses = string_to_process.split('#SPLIT#')
        c=0
        main_res=[]
        parentheses_res=[]
        while split_parentheses:
            sp=split_parentheses.pop(0).strip()
            if c%2==0:
                if sp: main_res.append(sp)
                c+=1
            else:
                if sp:parentheses_res.append([sp])
                if c<=n_parentheses: c+=1
        res=[main_res]
        res.extend(parentheses_res)
        res=[' '.join(i) for i in res]
        return res



    def final_processing(self, string_to_process):
        new_str = string_to_process.strip(' ')
        new_str=new_str.replace('[]','')
        new_str=new_str.replace('()','')
        new_str=new_str.replace('_','')
        new_str=new_str.replace('-like ',' ')
        new_str=new_str.replace('-domain-containing ',' ')
        new_str=new_str.replace('-associated ',' ')
        if 'biosynth' in new_str:new_str=new_str.replace('biosynth','synth')
        if new_str.endswith('-like'):new_str=new_str.replace('-like','')
        if new_str.endswith('-domain-containing'):new_str=new_str.replace('-domain-containing','')
        if new_str.endswith('-associated'):new_str=new_str.replace('-associated','')
        new_str = new_str.replace('. ', '!NEWLINE!')
        new_str = new_str.replace(' / ', '!NEWLINE!')
        new_str = new_str.replace('\t', '!NEWLINE!')
        new_str=self.split_by_alternatives(new_str)
        lines = new_str.split('!NEWLINE!')
        res = []
        for line in lines:
            #split_parentheses = line.split('#SPLIT#')
            #split_parentheses = [i for i in split_parentheses if i]
            split_parentheses=self.split_by_parentheses(line)
            for current_str in split_parentheses:
                current_str=self.remove_extra_parentheses(current_str)
                current_str = self.unite_terms(current_str)
                current_str = [i.strip('.') for i in current_str]
                current_str = [self.remove_extra_parentheses(i) for i in current_str if i]
                res.append(current_str)
        res=[[j.lower() for j in i if j] for i in res if i]
        return res


    def remove_go_obo_identifiers(self, sentence):
        '''
        ids to keep
        vz =  https://viralzone.expasy.org/
        SO = http://www.sequenceontology.org/
        metacyc
        reactome
        hgnc
        pfam
        chebi
        brenda

        cant keep KEGG because GO.obo has no distinction between KEGG's identifiers types...
        '''
        go_obo_pattern = '\[(goc|PR|CL|Wikipedia|CORUM|MetaCyc|GOC|ISBN|PMID|Reactome|CHEBI|GO|VZ|vz|gOC|HGNC|KEGG|KEGG_REACTION|UBERON|Pfam|RESID|MA|SO|UniPathway|MP|BRENDA|DOI|pmid|' \
                         'Wikipeda|Wikilpedia|MGI|DDANAT|PO|ABA):[A-Za-z\d\-]+(,\s(goc|PR|CL|Wikipedia|CORUM|MetaCyc|GOC|ISBN|PMID|Reactome|CHEBI|GO|VZ|vz|gOC|HGNC|KEGG|KEGG_REACTION|' \
                         'UBERON|Pfam|RESID|MA|SO|UniPathway|MP|BRENDA|DOI|pmid|Wikipeda|Wikilpedia|MGI|DDANAT|PO|ABA):[A-Za-z\d\-]+)*\]'
        res, go_obo_ids = self.remove_pattern(sentence, go_obo_pattern)
        go_obo_ids_res = {}
        http_pattern = re.compile('\[http.*\]')
        search_http = re.search(http_pattern, res)
        if search_http:
            res = res.replace(search_http.group(), '')
        if go_obo_ids:
            go_obo_ids = go_obo_ids.pop()
            go_obo_ids = go_obo_ids.replace('[', '')
            go_obo_ids = go_obo_ids.replace(']', '')
            go_obo_ids = go_obo_ids.split(', ')
            go_obo_ids_res = {}
            for i in go_obo_ids:
                go_obo_type, go_obo_id = i.split(':')
                if go_obo_type in ['vz', 'VZ', 'SO', 'MetaCyc', 'Reactome', 'HGNC', 'Pfam', 'CHEBI', 'BRENDA']:
                    if go_obo_type == 'vz':
                        go_obo_type = 'viralzone'
                    elif go_obo_type == 'VZ':
                        go_obo_type = 'viralzone'
                    elif go_obo_type == 'SO':
                        go_obo_type = 'seq_onto'
                    elif go_obo_type == 'BRENDA':
                        go_obo_type = 'enzyme_ec'
                    else:
                        go_obo_type = go_obo_type.lower()
                    if go_obo_type not in go_obo_ids_res:
                        go_obo_ids_res[go_obo_type] = set()
                    go_obo_ids_res[go_obo_type].add(go_obo_id)
        return res, go_obo_ids_res


    def remove_common_identifiers(self, sentence):
        res = ' ' + str(sentence) + ' '
        dict_ids = {'enzyme_ec': set(),
                    'tcdb': set(),
                    'kegg_ko': set(),
                    'tigrfam': set(),
                    'pfam': set(),
                    'cog': set(),
                    'go': set(),
                    'others': set()
                    }
        res, removed_ecs = self.remove_ecs(res, required_level=1)
        res = res.replace(' enzyme_ec:', '')
        dict_ids['enzyme_ec'].update(removed_ecs)
        tcdb_pattern = 'TC\s\d\.[A-Z\-](\.(\d+|\-)){1,2}'
        ko_pattern = 'K\d{4,}'
        tigrfam_pattern = 'TIGR\d{3,}'
        duf_pattern = '(DUF|duf)\d{2,}'
        pfam_pattern = '((U|u)?PF|pf)\d{3,}'
        cog_pattern = '(COG|cog)\d{3,}'
        go_pattern = 'GO:?\d{3,}'
        pubmed_pattern = 'PubMed:\d+'

        res, removed_ids = self.remove_pattern(res, tcdb_pattern)
        res = res.replace(' tcdb:', '')
        dict_ids['tcdb'].update(removed_ids)
        res, removed_ids = self.remove_pattern(res, tigrfam_pattern)
        res = res.replace(' tigrfam:', '')
        dict_ids['tigrfam'].update(removed_ids)
        res, removed_ids = self.remove_pattern(res, ko_pattern)
        res = res.replace(' kegg_ko:', '')
        dict_ids['kegg_ko'].update(removed_ids)
        res, removed_ids = self.remove_pattern(res, duf_pattern)
        res = res.replace(' pfam:', '')
        dict_ids['pfam'].update(removed_ids)
        res, removed_ids = self.remove_pattern(res, pfam_pattern)
        dict_ids['pfam'].update(removed_ids)
        res, removed_ids = self.remove_pattern(res, cog_pattern)
        res = res.replace(' cog:', '')
        dict_ids['cog'].update(removed_ids)
        res, removed_ids = self.remove_pattern(res, go_pattern)
        res = res.replace(' go:', '')
        removed_ids = {re.search('\d+', i).group() for i in removed_ids}
        dict_ids['go'].update(removed_ids)
        res, removed_ids = self.remove_pattern(res, pubmed_pattern)
        res = res.replace(' pubmed:', '')
        res, go_obo_identifiers = self.remove_go_obo_identifiers(res)
        for go_obo_type in go_obo_identifiers:
            if go_obo_type == 'enzyme_ec':
                dict_ids['enzyme_ec'].update(go_obo_identifiers[go_obo_type])
            else:
                dict_ids[go_obo_type] = go_obo_identifiers[go_obo_type]
        return res, dict_ids



    def pre_process_string(self, sentence):
        if not sentence: return [], []
        res,dict_ids = self.remove_common_identifiers(sentence)
        res = self.replace_punctuation(res)

        digit_pattern = '\d+(\.\d+)?'
        res, _ = self.remove_pattern(res, digit_pattern)
        for bio_pattern in [
            '[A-Z\d]+_[A-Z\d]+(_[A-Z\d]+)?',
            '[A-Z]+\d{3,}(\.\d+)?([A-Z]+)?',
            '[a-z]{3,}([A-Z]|\d+)(\d+)?',
            '[A-Z]{3,}(\d+)?',
            '[A-Z]+[a-z]+[A-Z]+',
            '[A-Z]\d+-\d+',
            '[A-Z][a-z]+\d+',
            '\d[A-Z]',
        ]:
            res, ids_removed = self.remove_pattern(res, bio_pattern)
            dict_ids['others'].update(ids_removed)


        # cleaning blank parenthesis
        res, _ = self.remove_pattern(res, '\(\s+\)')
        res, _ = self.remove_pattern(res, '\[\s+\]')
        res = self.remove_bad_pattern_ions(res)
        res = self.simplify_ions(res)
        res = self.replace_roman_numerals(res)
        res = self.remove_parentheses(res)
        all_ids = set()
        for id_type in dict_ids:
            if id_type == 'others':
                all_ids.update(dict_ids['others'])
            else:
                for id_str in dict_ids[id_type]:
                    all_ids.add(id_type + ':' + id_str)
        proper_tokens=self.final_processing(res)
        return proper_tokens, all_ids


class WordNetTagger(SequentialBackoffTagger):
    def __init__(self,perceptron_tagger,go_terms=None, *args, **kwargs):
        SequentialBackoffTagger.__init__(self, *args, **kwargs)
        self.stemmer = PorterStemmer()
        #for universal tagger
        self.placeholder = 'XXXXX'
        self.perceptron_tagger=perceptron_tagger
        self.go_terms=set()
        for g in go_terms:
            #these tags were manually reviewed so they wouldnt intefere with go terms tagging
            if self.tag_tokens_perceptron([g]) not in ['ADP','CONJ','DET','PRON','PRT']:
                self.go_terms.add(g)
        self.wordnet_tag_map = {'NOUN': 'NOUN','ADJ': 'ADJ',
                                'ADV': 'ADV','VERB': 'VERB',
                                #placeholder for entities
                                self.placeholder:self.placeholder}

    def tag_tokens_perceptron(self, tokens):
        return _pos_tag(tokens, tagset='universal', tagger=self.perceptron_tagger,lang='eng')

    def choose_tag(self, tokens, index, history):
        word = tokens[index]
        word=word.strip()
        if word ==self.placeholder: return self.placeholder
        #adding terms from pre-processing
        if len(word.split(' '))>1:
            terms=word.split(' ')
            word=terms[0]
        fd = FreqDist()
        for synset in wn.synsets(word):
            lex_name=synset.lexname().split('.')[0]
            fd[lex_name] += 1
        wordnet_res,go_res=None,None
        if fd.keys(): wordnet_res= self.wordnet_tag_map.get(fd.max())
        tagger_res=self.tag_tokens_perceptron([word])[0][1]
        if word in self.go_terms: go_res= 'NOUN'
        if wordnet_res: return wordnet_res
        elif tagger_res:
            if tagger_res!='NOUN':
                return tagger_res
        elif go_res: return go_res
        else: return None

class Metadata():

    # this unzips to the same directory!
    def gunzip(self,source_filepath, dest_filepath=None, block_size=65536, remove_source=True):
        if not dest_filepath:
            dest_filepath = source_filepath.strip('.gz')
        if os.path.isdir(dest_filepath):
            file_name = source_filepath.split(SPLITTER)[-1].replace('.gz', '')
            dest_filepath = add_slash(dest_filepath) + file_name
        print('Gunzipping ', source_filepath, 'to', dest_filepath, flush=True)
        with gzip_open(source_filepath, 'rb') as s_file, \
                open(dest_filepath, 'wb') as d_file:
            while True:
                block = s_file.read(block_size)
                if not block:
                    break
                else:
                    d_file.write(block)
            d_file.write(block)
        if remove_source: os.remove(source_filepath)

    def download_file_http(self, url, file_path, c):
        if c > 5:
            self.download_file_http_failsafe(url, file_path)
        else:
            with requests.get(url, stream=True) as r:
                with open(file_path, 'wb') as f:
                    shutil.copyfileobj(r.raw, f)

    # slower but safer
    def download_file_http_failsafe(self, url, file_path):
        with requests.Session() as session:
            get = session.get(url, stream=True)
            if get.status_code == 200:
                with open(file_path, 'wb') as f:
                    for chunk in get.iter_content(chunk_size=1024):
                        f.write(chunk)

    def download_file_ftp(self, url, file_path):
        with closing(request.urlopen(url)) as r:
            with open(file_path, 'wb') as f:
                shutil.copyfileobj(r, f)

    def download_file(self, url, output_folder='', retry_limit=10):
        file_path = output_folder + url.split('/')[-1]
        try:
            target_file = request.urlopen(url)
        except:
            print('Cannot download target url', url)
            return
        target_size = target_file.info()['Content-Length']
        transfer_encoding = target_file.info()['Transfer-Encoding']
        if target_size: target_size = int(target_size)
        if os.path.exists(file_path):
            if transfer_encoding == 'chunked':
                return
            elif os.stat(file_path).st_size == target_size:
                print('Not downloading from ' + url + ' since file was already found!', flush=True)
                return
            else:
                os.remove(file_path)
        print('Downloading from ' + url + '. The file will be kept in ' + output_folder, flush=True)
        c = 0
        while c <= retry_limit:
            if 'ftp' in url:
                self.download_file_ftp(url, file_path)
            else:
                self.download_file_http(url, file_path, c)
            if transfer_encoding == 'chunked': return
            if os.stat(file_path).st_size == target_size: return
            c += 1
        print('Did not manage to download the following url correctly:\n' + url)
        raise Exception

    def get_taxon_ids_NOGT(self, url):
        webpage = None
        c = 0
        while not webpage and c <= 10:
            req = requests.get(url)
            try:
                webpage = req.text
            except:
                c += 1
        taxons_search = re.findall('href="\d+/"', webpage)
        taxons = [re.search('\d+', i).group() for i in taxons_search]
        return taxons

    def download_all_metadata(self):
        pfam_metadata = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz'
        ko_list = 'ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz'
        ncbi_metadata = 'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv'


        for url in [pfam_metadata,ko_list,ncbi_metadata]:
            self.download_file(url,output_folder=self.metadata_folder)
        self.download_file('http://purl.obolibrary.org/obo/go.obo',output_folder=self.unifunc_resources_folder)

        eggnog_downloads_page = 'http://eggnog5.embl.de/download/latest/per_tax_level/'
        for taxon_id in self.get_taxon_ids_NOGT(eggnog_downloads_page):
            url = 'http://eggnog5.embl.de/download/latest/per_tax_level/' + str(taxon_id) + '/' + str(taxon_id) + '_annotations.tsv.gz'
            self.download_file(url,output_folder=self.metadata_folder)
        for file in os.listdir(self.metadata_folder):
            if file.endswith('.gz'):
                self.gunzip(source_filepath=self.metadata_folder+file,dest_filepath=self.metadata_folder+file.strip('.gz'))

    def build_freq_dict_uniprot_reference(self):
        file_reference=self.unifunc_resources_folder + 'uniprot.tab'
        print('Adding terms to frequency counter from ', file_reference, flush=True)
        stop_words = set(stopwords.words("english"))
        with open(file_reference) as file:
            line = file.readline()
            line = file.readline()
            while line:
                line = line.strip('\n')
                func_res = self.process_uniprot(line)
                if func_res:
                    self.document_counter += 1
                    query, annotation, function_description = func_res
                    for desc in [annotation, function_description]:
                        already_added = set()
                        # each annotation line is a document,each word a token, the whole collection of annotations is the corpus
                        list_words = self.pre_process_string(desc)[0]
                        for sw in list_words:
                            n_grams = self.generate_n_grams(sw)
                            # for word in sw:
                            for word in n_grams:
                                lower_word = word.lower()
                                if self.is_good_word(lower_word, stop_words):
                                    if lower_word not in self.word_counter: self.word_counter[lower_word] = 0
                                    if lower_word not in already_added:
                                        self.word_counter[lower_word] += 1
                                        already_added.add(lower_word)
                line = file.readline()

    def add_to_word_counter(self,split_tokens):
        stop_words = set(stopwords.words("english"))
        if split_tokens:
            already_added = set()
            self.document_counter += 1
            for sw in split_tokens:
                n_grams = self.generate_n_grams(sw)
                # for word in sw:
                for word in n_grams:
                    lower_word = word.lower()
                    if self.is_good_word(lower_word, stop_words):
                        if lower_word not in self.word_counter: self.word_counter[lower_word] = 0
                        if lower_word not in already_added:
                            self.word_counter[lower_word] += 1
                            already_added.add(lower_word)

    def build_freq_dict_go_terms(self):
        print('Adding terms to frequency counter from ', self.go_terms_path, flush=True)
        with open(self.go_terms_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                split_tokens = None
                if line:
                    if 'name: ' in line[0:6] and 'obsolote' not in line.lower():
                        go_name = line.replace('name: ', '')
                        split_tokens = self.pre_process_string(go_name)[0]
                    elif 'synonym: ' in line[0:9] and 'EXACT' in line and 'obsolote' not in line.lower():
                        go_syn = line.replace('synonym: ', '')
                        go_syn = go_syn.replace('EXACT ', '')
                        go_syn = go_syn.replace('[]', '')
                        split_tokens = self.pre_process_string(go_syn)[0]
                    elif 'def: ' in line[
                                    0:5] and 'Catalysis of the reaction' not in line and 'obsolote' not in line.lower():
                        go_def = line.replace('def: ', '')
                        go_def = go_def.split('\"')[1]
                        split_tokens = self.pre_process_string(go_def)[0]
                    self.add_to_word_counter(split_tokens)
                line = file.readline()

    def build_freq_dict_pfam(self):
        file_reference=self.metadata_folder+'Pfam-A.hmm.dat'
        print('Adding terms to frequency counter from ', file_reference, flush=True)
        with open(file_reference) as file:
            line = file.readline()
            while line:
                line = line.strip('\n').split('   ')
                if len(line) == 2:
                    row_header, row_description = line
                    if row_header == '#=GF ID':
                        stop = True

                    elif row_header == '#=GF DE' and stop:
                        stop = False
                        split_tokens = self.pre_process_string(row_description)[0]
                        self.add_to_word_counter(split_tokens)
                line = file.readline()

    def build_freq_dict_kofam(self):
        file_reference=self.metadata_folder+'ko_list'
        print('Adding terms to frequency counter from ', file_reference, flush=True)
        with open(file_reference) as file:
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                description = line[-1]
                if '[EC:' in description:
                    description, temp_links = description.split('[EC:')
                split_tokens = self.pre_process_string(description)[0]
                self.add_to_word_counter(split_tokens)
                line = file.readline()

    def build_freq_dict_pfm(self):
        file_reference = self.metadata_folder + 'hmm_PGAP.tsv'
        print('Adding terms to frequency counter from ', file_reference, flush=True)

        with open(file_reference) as file:
            line = file.readline()
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                description=line[10]
                description = description.replace('(Provisional)', '')
                split_tokens = self.pre_process_string(description)[0]
                self.add_to_word_counter(split_tokens)
                line = file.readline()

    def build_freq_dict_eggnog(self):
        for f in os.listdir(self.metadata_folder):
            if f.endswith('_annotations.tsv'):
                file_reference = self.metadata_folder + f
                print('Adding terms to frequency counter from ', file_reference, flush=True)
                with open(file_reference) as file:
                    line = file.readline()
                    while line:
                        line = line.strip('\n').split('\t')
                        description = line[3]
                        split_tokens = self.pre_process_string(description)[0]
                        self.add_to_word_counter(split_tokens)
                        line = file.readline()

    def build_freq_dict_all(self):
        self.build_freq_dict_uniprot_reference()
        self.build_freq_dict_go_terms()
        self.build_freq_dict_pfam()
        self.build_freq_dict_kofam()
        self.build_freq_dict_pfm()
        self.build_freq_dict_eggnog()










class UniFunc(Pre_Processer, Word_Weighter, Metadata):
    def __init__(self,compile_metadata=False):
        #this is just used to trigger gene ontologies look up, might remove it
        self.nlp_threshold = 0.9
        self.n_grams_range = [1]
        str_n_gram = '_'.join([str(i) for i in self.n_grams_range])
        self.unifunc_resources_folder=unifunc_folder+'Resources/'
        self.go_terms_path=self.unifunc_resources_folder+'go.obo'
        self.verbs_txt=self.unifunc_resources_folder+'verbs.tab'
        self.word_counter_pickled = self.unifunc_resources_folder+'frequency_dict_n_grams_'+str_n_gram+'.pickle'
        self.verbs_set_pickled = self.unifunc_resources_folder+'verbs_set.pickle'
        if compile_metadata:
            self.metadata_folder=self.unifunc_resources_folder+'Metadata_dbs/'
            if not os.path.exists(self.metadata_folder) and not os.path.exists(self.word_counter_pickled):
                Path(self.metadata_folder).mkdir(parents=True, exist_ok=True)
                self.download_all_metadata()

        self.download_nltk_resources()
        #uses penn treebank corpus
        self.tagger = PerceptronTagger()

        self.word_counter={}
        self.verbs_set={}
        self.document_counter=0
        self.good_identifiers={'enzyme_ec','tcdb','kegg_ko','tigrfam','pfam','cog','go','viralzone','seq_onto'}
        self.words_to_remove = ['mainrole','sub1role','protein','proteins',
                                'enzyme','enzymes','putative','activity',
                                'process','unknown','function','functions',
                                'processes','responsible','probable',
                                'other','complex','integral',
                                'metabolic','identical','type',
                                'related','large','small',
                                'accessory','major','related,'
                                'variable','potential','specific',
                                'hypothetical','particle','particles',
                                ]
        self.build_frequency_dict()
        self.get_verbs_set()
        self.pickled_go_syns = self.go_terms_path+'.pickle_syns'
        self.pickled_go_terms = self.go_terms_path+'.pickle_terms'
        self.pickled_go_dict = self.go_terms_path+'.pickle_dict'
        self.go_syns=set()
        self.go_terms=set()
        self.go_dict=dict()
        self.parse_go_terms()
        self.wordnet_tagger = WordNetTagger(go_terms=self.go_terms, perceptron_tagger=self.tagger)
        self.tags = {}
        self.default_word_count=median([self.word_counter[i] for i in self.word_counter if self.word_counter[i]>1])


    def tag_tokens_perceptron(self, tokens):
        return _pos_tag(tokens, tagset='universal', tagger=self.tagger, lang='eng')

    def tag_tokens_wordnet(self, tokens):
        return self.wordnet_tagger.tag(tokens)

    def __str__(self):
            res='Tags:\n'
            for t in self.tags:
                res+=t+': '+str(self.tags[t])+'\n'
            return res

    def download_nltk_resources(self):
        try:
            nltk_download('stopwords',quiet=True)
        except:
            print('Already downloaded stopwords')
        try:
            nltk_download('averaged_perceptron_tagger',quiet=True)
        except:
            print('Already downloaded Perceptron tagger')
        try:
            nltk_download('universal_tagset',quiet=True)
        except:
            print('Already downloaded Universal tagset!')
        try:
            nltk_download('wordnet',quiet=True)
        except:
            print('Already downloaded Wordnet!')

    ####GO SCORING
    def save_go_pickle(self):
        with open(self.pickled_go_syns, 'wb') as handle:
            pickle.dump(self.go_syns, handle,protocol=4)
        with open(self.pickled_go_terms, 'wb') as handle:
            pickle.dump(self.go_terms, handle,protocol=4)
        with open(self.pickled_go_dict, 'wb') as handle:
            pickle.dump(self.go_dict, handle,protocol=4)



    def load_go_pickle(self):
        if os.path.exists(self.pickled_go_syns):
            with open(self.pickled_go_syns, 'rb') as handle:
                self.go_syns = pickle.load(handle)
        if os.path.exists(self.pickled_go_terms):
            with open(self.pickled_go_terms, 'rb') as handle:
                self.go_terms = pickle.load(handle)
        if os.path.exists(self.pickled_go_dict):
            with open(self.pickled_go_dict, 'rb') as handle:
                self.go_dict = pickle.load(handle)


    def token_list_too_small(self,tokens_list,all_tokens_list):
        res=[]
        if len(tokens_list)>1: return False
        for tl in all_tokens_list:
            res.append(len(tl))
        max_len= max(res)
        if len(tokens_list)==1 and len(tokens_list)<max_len:
            return True
        return False

    def parse_go_terms(self):
        self.load_go_pickle()
        if not self.go_syns or not self.go_terms or not self.go_dict:
            print('Parsing GO terms with ',self.go_terms_path,flush=True)
            go_terms=set()
            go_dict={}
            with open(self.go_terms_path) as file:
                line=file.readline()
                while line:
                    line=line.strip('\n')
                    if line:
                        if 'id: ' in line[0:4]:
                            go_id=line.replace('id: ','')
                            go_id=go_id.lower()
                            go_dict[go_id]={'synonyms':set(),'identifiers':set()}
                        elif 'alt_id: ' in line[0:8]:
                            alt_go_id=line.replace('alt_id: ','')
                            alt_go_id=alt_go_id.lower()
                            go_dict[alt_go_id]=go_dict[go_id]


                        elif 'name: ' in line[0:6]:
                            go_name=line.replace('name: ','')
                            split_tokens,all_ids=self.pre_process_string(go_name)
                            for i in split_tokens:
                                tags=self.tag_tokens_perceptron(i)
                                for t in tags:
                                    if t[1] not in ['ADP','CONJ','DET','PRON','PRT'] and\
                                        t[0].lower() not in self.words_to_remove and\
                                        len(t[0])>1 and\
                                        not  re.search('\d+',t[0]):
                                        go_terms.add(t[0])
                            for token_list in split_tokens:
                                if not self.token_list_too_small(token_list,split_tokens):
                                    go_dict[go_id]['synonyms'].add(' '.join(token_list))
                            go_dict[go_id]['identifiers'].update(all_ids)
                        elif 'synonym: ' in line[0:9] and 'EXACT' in line:
                            go_syn=line.replace('synonym: ','')
                            go_syn=go_syn.replace('EXACT ','')
                            go_syn=go_syn.replace('[]','')

                            split_tokens,all_ids=self.pre_process_string(go_syn)
                            for token_list in split_tokens:
                                if not self.token_list_too_small(token_list,split_tokens):
                                    go_dict[go_id]['synonyms'].add(' '.join(token_list))
                            go_dict[go_id]['identifiers'].update(all_ids)

                    line = file.readline()
                    if '[Typedef]' in line: line=None
            for go_id in go_dict:
                self.go_syns.add(frozenset(go_dict[go_id]['synonyms']))
            self.go_terms=go_terms
            self.go_dict=go_dict
            self.save_go_pickle()

    def has_go_match(self,test_syn,ref_syn):
        for syn_set in self.go_syns:
            if ref_syn in syn_set and test_syn in syn_set:
                return True
        return False



    ####NLP SCORING


    def remove_unwanted_tokens(self,tagged_tokens):
        '''
        POS tag list:
        for the universal tagset:
            ADJ 	adjective 	new, good, high, special, big, local
            ADP 	adposition 	on, of, at, with, by, into, under
            ADV 	adverb 	really, already, still, early, now
            CONJ 	conjunction 	and, or, but, if, while, although
            DET 	determiner, article 	the, a, some, most, every, no, which
            NOUN 	noun 	year, home, costs, time, Africa
            NUM 	numeral 	twenty-four, fourth, 1991, 14:24
            PRT 	particle 	at, on, out, over per, that, up, with
            PRON 	pronoun 	he, their, her, its, my, I, us
            VERB 	verb 	is, say, told, given, playing, would
            . 	punctuation marks 	. , ; !
            X 	other 	ersatz, esprit, dunno, gr8, univeristy
        '''
        res = []
        #ideally we'd only keep the nouns and numbers (maybe verbs?)but there's a lot of false positives...
        tags_to_remove=[self.wordnet_tagger.placeholder,'DET','PRON','PRT','CONJ','ADP']
        #from tigrfam
        stop_words = set(stopwords.words("english"))
        for t in tagged_tokens:
            if t[1] not in ['NOUN']:
                if t[1] not in self.tags: self.tags[t[1]]=set()
                self.tags[t[1]].add(t[0])
            if t[1] not in tags_to_remove and len(t[0])>1 and t[0].lower() not in self.words_to_remove and t[0] not in stop_words and t[0] not in self.verbs_set:
                res.append(t[0])

        return res

    def choose_best_tagging(self,wornet_tagging,default_tagging):
        res=[]
        for word_tag_i in range(len(wornet_tagging)):
            if default_tagging[word_tag_i][1]=='NUM':res.append(default_tagging[word_tag_i])
            elif wornet_tagging[word_tag_i][1]: res.append(wornet_tagging[word_tag_i])
            else: res.append(default_tagging[word_tag_i])
        return res




    def process_string_nlp(self,tokens_lists):
        res=[]
        for tokens in tokens_lists:
            # word lemmatization - doesnt make much sense since we are not trying to classify text meaning. It won't change the tokens that much and might even produce errors
            # word stemming - makes a bit more sense but once again it also produces weird errors
            temp_tokens=[i for i in tokens if i]
            default_tagged_tokens = self.tag_tokens_perceptron(temp_tokens)
            #_ = self.tag_tokens_wordnet(temp_tokens)
            wordnet_tagged_tokens = self.tag_tokens_wordnet(temp_tokens)
            tagged_tokens=self.choose_best_tagging(wordnet_tagged_tokens,default_tagged_tokens)
            removed_tags = self.remove_unwanted_tokens(tagged_tokens)
            if removed_tags:
                res.append(removed_tags)
        return res


    def get_best_syn(self,original_syns_set,original_list_words,edit_dist_perc=0.05):
        list_words=set([i.lower() for i in original_list_words])
        syns_set=set([i.lower() for i in original_syns_set])
        for syn in syns_set:
            if syn in list_words: return syn
        for syn in syns_set:
            for word in list_words:
                max_len=max([len(syn),len(word)])
                if edit_distance(syn,word)/max_len<=edit_dist_perc:
                    return word

    def add_wordnet_lemmas(self,res,w):
        for synset in wn.synsets(w):
            for lemma in synset.lemmas():
                to_add = lemma.name()
                if '_' not in to_add:
                    res.add(lemma.name())

    def add_alternative_tokens(self,vector1):
        #this will mostly add alternative suffixes
        res={}
        for w in vector1:
            res[w]={w}
            self.add_wordnet_lemmas(res[w],w)
            if w=='decarboxylase': res[w].add('carboxy-lyase')
            elif w=='carboxy-lyase': res[w].add('decarboxylase')
            #latin plurals
            if w.endswith('ae'):        res[w].add(w[:-1])
            if w.endswith('er'):        res[w].add(w[:-2])
            if w.endswith('exes'):      res[w].add(w[:-2])
            if w.endswith('ices'):      res[w].add(w[:-4]+'ex')
            if w.endswith('eaus'):      res[w].add(w[:-1])
            if w.endswith('eaux'):      res[w].add(w[:-1])
            if w.endswith('ia'):        res[w].add(w[:-1]+'on')
            if w.endswith('ions'):      res[w].add(w[:-1])
            if w.endswith('es'):        res[w].add(w[:-2]+'is')
            if w.endswith('os'):        res[w].add(w[:-1])
            if w.endswith('oes'):       res[w].add(w[:-2])
            if w.endswith('uses'):      res[w].add(w[:-2])
            if w.endswith('i'):         res[w].add(w[:-1]+'us')
            #normal plurals
            if w.endswith('s'):         res[w].add(w[:-1])
            #other suffixes
            if w.endswith('ation'):     res[w].add(w[:-5 ]+'e')
            if w.endswith('ation'):     res[w].add(w[:-3 ]+'e')
            if w.endswith('esis'):      res[w].add(w[:-3 ]+'tic')
            if w.endswith('etic'):      res[w].add(w[:-3 ]+'sis')
            #other latin plurals
            if w.endswith('a'):         res[w].add(w[:-1]+'um')
            if w.endswith('a'):         res[w].add(w[:-1]+'on')
            if w.endswith('ic'):        res[w].add(w[:-1]+'on')
            if w.endswith('i'):         res[w].add(w[:-1]+'on')
            #remove parentheses
            if w.startswith('[') and w.endswith(']'): res[w].add(w.strip('[]'))
            if w.startswith('(') and w.endswith(')'): res[w].add(w.strip('()'))
        return res

    def find_matching_tokens(self,synset1,synset2):
        #all words
        matches_found={i:None for i in set(synset1.keys()).union(synset2.keys())}
        for w1 in synset1:
            if not matches_found[w1]:
                synsw1=synset1[w1]
                for w2 in synset2:
                    if not matches_found[w2]:
                        if w1==w2:
                            matches_found[w1]=w2
                            matches_found[w2]=w1
                        else:
                            synsw2 = synset2[w2]
                            syns_intersect=synsw1.intersection(synsw2)
                            if syns_intersect:
                                best_syn=syns_intersect.pop()
                                matches_found[w1]=best_syn
                                matches_found[w2]=best_syn
        return matches_found

    def replace_tokens_with_alternatives(self,tokens_list,matches_found):
        res=[]
        for i in range(len(tokens_list)):
            if matches_found[tokens_list[i]]: res.append(matches_found[tokens_list[i]])
            else: res.append(tokens_list[i])
        return res

    def get_matching_suffixes(self,vector1,vector2):
        vector1=[i.lower() for i in vector1]
        vector2=[i.lower() for i in vector2]
        synset1=self.add_alternative_tokens(vector1)
        synset2=self.add_alternative_tokens(vector2)
        matches_found=self.find_matching_tokens(synset1,synset2)
        improved_vector1=self.replace_tokens_with_alternatives(vector1,matches_found)
        improved_vector2=self.replace_tokens_with_alternatives(vector2,matches_found)
        return improved_vector1,improved_vector2




    def build_vector(self,iterable1, iterable2):
        counter1 = self.calculate_tf_idf(iterable1)
        counter1= self.calculate_scaled_tf_idf(counter1)
        counter2 = self.calculate_tf_idf(iterable2)
        counter2= self.calculate_scaled_tf_idf(counter2)
        all_items = set(counter1.keys()).union(set(counter2.keys()))
        vector1,vector2=[],[]
        for word in all_items:
            if word in counter1:vector1.append(counter1[word])
            else: vector1.append(0)
            if word in counter2:vector2.append(counter2[word])
            else: vector2.append(0)
        return vector1, vector2

    def jaccard_distance(self,label1, label2):
        union_labels=len(label1.union(label2))
        intersection_labels=len(label1.intersection(label2))
        if not union_labels: return 1
        return (union_labels - intersection_labels) / union_labels

    def score_text(self,vector_1,vector_2):
        if not vector_1 or not vector_2: return 0
        if self.wordnet_tagger.placeholder.lower() in vector_1: vector_1.remove(self.wordnet_tagger.placeholder.lower())
        if self.wordnet_tagger.placeholder.lower() in vector_2: vector_2.remove(self.wordnet_tagger.placeholder.lower())
        #we do this so late because
        improved_vector_1,improved_vector_2 = self.get_matching_suffixes(vector_1,vector_2)
        v1,v2=self.build_vector(improved_vector_1,improved_vector_2)
        if not any(v1) or not any(v2):            text_score=0
        else:            text_score=1-cosine_distance(v1,v2)
        return text_score

    def check_token_intersection(self,tokens_lists_1,tokens_lists_2):
        for t1 in tokens_lists_1:
            for t2 in tokens_lists_2:
                if set(t1).intersection(set(t2)): return True
        return False

    def get_output(self,score,redirect_verbose,console_output,threshold,only_return=False):
        if only_return: return round(score,5)
        if threshold is not None:
            if score>=threshold:
                if console_output:
                    print(str(True) , flush=True, file=redirect_verbose)
                    redirect_verbose.close()
                print(str(True))
                return True
            else:
                if console_output:
                    print(str(False), flush=True, file=redirect_verbose)
                    redirect_verbose.close()
                print(str(False))
                return False
        if console_output:
            print(str(round(score,5)), flush=True,file=redirect_verbose)
            redirect_verbose.close()
        print(str(round(score, 5)))
        return round(score,5)


    def extract_bad_identifier_from_tokens(self,bad_identifiers,tokens_list):
        all_tokens = set()
        for tl in tokens_list:
            for token in tl: all_tokens.add(token.lower())
        res=[]
        temp_bad_identifiers=set([i.lower() for i in bad_identifiers])
        to_remove=temp_bad_identifiers.intersection(all_tokens)
        for tl in tokens_list:
            temp=[]
            for token in tl:
                if token.lower() not in to_remove: temp.append(token)
            res.append(temp)
        return res,to_remove

    def clean_tokens_list_of_identifiers(self,bad_identifiers_1,bad_identifiers_2,tokens_lists_1,tokens_lists_2):
        temp_tokens_list_1,temp_bad_identifiers_1 = self.extract_bad_identifier_from_tokens(bad_identifiers_2,tokens_lists_1)
        temp_bad_identifiers_1.update(bad_identifiers_1)
        temp_bad_identifiers_1=set([i.upper() for i in temp_bad_identifiers_1])
        temp_tokens_list_2,temp_bad_identifiers_2 = self.extract_bad_identifier_from_tokens(bad_identifiers_1,tokens_lists_2)
        temp_bad_identifiers_2.update(bad_identifiers_2)
        temp_bad_identifiers_2=set([i.upper() for i in temp_bad_identifiers_2])
        return temp_tokens_list_1,temp_tokens_list_2,temp_bad_identifiers_1,temp_bad_identifiers_2

    def get_similarity_score(self,string_1, string_2,console_output=None,verbose=False,only_text=False,threshold=None,only_return=False):
        if console_output: redirect_verbose=open(console_output,'a+')
        else:redirect_verbose=None

        if verbose:
            print('Analysing strings:',flush=True,file=redirect_verbose)
            print('\tString 1: ',string_1,flush=True,file=redirect_verbose)
            print('\tString 2: ',string_2,flush=True,file=redirect_verbose)
        if not string_1 or not string_2: return self.get_output(score=0,redirect_verbose=redirect_verbose,console_output=console_output,threshold=threshold,only_return=only_return)
        temp_string_1=str(string_1)
        temp_string_2=str(string_2)
        #NLP SCORING
        #we dont currently take into account parentheses information (except if its an abbreviation), as it often adds counfounders
        tokens_lists_1,ids_removed_1 = self.pre_process_string(temp_string_1)
        tokens_lists_2,ids_removed_2 = self.pre_process_string(temp_string_2)
        if verbose:
            print('We pre-processed our strings, here are the tokens and identifiers for each string:',flush=True,file=redirect_verbose)
            print('\tList of tokens for string 1:',tokens_lists_1,flush=True,file=redirect_verbose)
            print('\tList of tokens for string 2:',tokens_lists_2,flush=True,file=redirect_verbose)
            print('\tList of identifiers for string 1:',ids_removed_1,flush=True,file=redirect_verbose)
            print('\tList of identifiers for string 2:',ids_removed_2,flush=True,file=redirect_verbose)

        #for identifiers we know are highly resolved
        good_identifiers_1={i for i in ids_removed_1 if i.split(':')[0] in self.good_identifiers}
        good_identifiers_2={i for i in ids_removed_2 if i.split(':')[0] in self.good_identifiers}
        bad_identifiers_1={i for i in ids_removed_1 if i not in good_identifiers_1}
        bad_identifiers_2={i for i in ids_removed_2 if i not in good_identifiers_2}
        if verbose:
            print("We found these common database identifiers, we will now check if they intersect.",flush=True,file=redirect_verbose)
            print('\tList of common database identifiers for string 1:', good_identifiers_1,flush=True,file=redirect_verbose)
            print('\tList of common database identifiers for string 2:', good_identifiers_2,flush=True,file=redirect_verbose)
        if not only_text:
            if good_identifiers_1.intersection(good_identifiers_2): return self.get_output(score=1,redirect_verbose=redirect_verbose,console_output=console_output,threshold=threshold,only_return=only_return)
        if verbose: print('There was no intersection between any of the common database identifiers, so we will keep analysing our strings!',flush=True,file=redirect_verbose)
        #to avoid using the taggers unnecessarily we first check if there are any intersections
        token_intersection = self.check_token_intersection(tokens_lists_1,tokens_lists_2)
        tokens_lists_1,tokens_lists_2,bad_identifiers_1,bad_identifiers_2=self.clean_tokens_list_of_identifiers(bad_identifiers_1,bad_identifiers_2,tokens_lists_1,tokens_lists_2)
        bad_ids_intersection = bad_identifiers_1.intersection(bad_identifiers_2)
        if not token_intersection and not bad_ids_intersection:
            if verbose:
                print('There was no intersection of tokens so we will stop here!',flush=True,file=redirect_verbose)
            return self.get_output(score=0,redirect_verbose=redirect_verbose,console_output=console_output,threshold=threshold,only_return=only_return)
        if verbose:   print('There was an intersection of tokens so we will keep analysing the strings!',flush=True,file=redirect_verbose)
        vector_1_list = self.process_string_nlp(tokens_lists_1)
        vector_2_list = self.process_string_nlp(tokens_lists_2)
        if verbose:
            print('We now tagged our tokens and removed unwanted tokens, here are the tokens for each string:',flush=True,file=redirect_verbose)
            print('\tList of tokens for string 1:',vector_1_list,flush=True,file=redirect_verbose)
            print('\tList of tokens for string 2:',vector_2_list,flush=True,file=redirect_verbose)
        text_score_list=[]
        min_len=min(len(vector_1_list),len(vector_2_list))
        all_vector_1=set()
        all_vector_2=set()
        for vector_1 in vector_1_list:
            all_vector_1.update(vector_1)
            for vector_2 in vector_2_list:
                all_vector_2.update(vector_2)
                current_score=self.score_text(vector_1,vector_2)
                text_score_list.append(current_score)
                if verbose:
                    print('Finding similarity between:',flush=True,file=redirect_verbose)
                    print('\t',vector_1,flush=True,file=redirect_verbose)
                    print('\t',vector_2,flush=True,file=redirect_verbose)
                    print('Similarity score between these two vectors was ',round(current_score,5),flush=True,file=redirect_verbose)
        #some annotations are split into multiple sentences when they shouldn't. so we merge all the tokens and check their similarity here
        if all_vector_1 and all_vector_2:
            if len(all_vector_1)==len(vector_1_list[0]) and len(vector_1_list)==1 and\
                    len(all_vector_2)==len(vector_2_list[0]) and len(vector_2_list)==1:
                pass
            else:
                all_vector_score = self.score_text(all_vector_1, all_vector_2)
                text_score_list.append(all_vector_score)
                if verbose:
                    print('Finding similarity between:', flush=True, file=redirect_verbose)
                    print('\t', all_vector_1, flush=True, file=redirect_verbose)
                    print('\t', all_vector_2, flush=True, file=redirect_verbose)
                    print('Similarity score between these two vectors was ', round(all_vector_score, 5),
                          flush=True,file=redirect_verbose)

        if text_score_list:
            wanted_top = min_len
            text_score_list=sorted(text_score_list)[-wanted_top:]
            text_score=mean(text_score_list)
        else: text_score=0
        ids_removed_score= 1-self.jaccard_distance(set(bad_identifiers_1),set(bad_identifiers_2))
        if verbose:
            if bad_identifiers_1 and bad_identifiers_2:
                print('Calculate the Jaccard distance between:',flush=True,file=redirect_verbose)
                print('\t',bad_identifiers_1,flush=True,file=redirect_verbose)
                print('\t',bad_identifiers_2,flush=True,file=redirect_verbose)
                print('The Jaccard distance between identifiers is:',round(ids_removed_score,5),flush=True,file=redirect_verbose)
            else: print('No IDs to compare so we skipped Jaccard distance calculation',flush=True,file=redirect_verbose)
        #more weight for the entities
        if not only_text:
            #if no tokens to score but has ids
            if ids_removed_score and not text_score_list:                       score=ids_removed_score
            #if score ids is sufficiently high, and tokens score is very high
            elif ids_removed_score>=0.5 and text_score==1:                      score=text_score
            #if ids score is above 0
            elif ids_removed_score: score=(ids_removed_score+text_score)/3
            elif bad_identifiers_1 and bad_identifiers_2 and not ids_removed_score: score=(ids_removed_score+text_score)/2
            else: score=text_score
        else:
            score=text_score
        if not only_text:
            #GO SCORING
            if score<self.nlp_threshold and score >0.3:
                if verbose:
                    print('The similarity score between our strings was too low so we will try to find a match between the strings and a gene ontology annotation.',flush=True,file=redirect_verbose)
                if self.has_go_match(string_1.lower(),string_2.lower()):
                    if verbose:
                        print('Found an intersection between our strings and the gene ontologies database. The strings are the same!',flush=True,file=redirect_verbose)
                    go_score=1
                else:
                    if verbose:
                        print('Didnt find any intersection between our strings and the gene ontologies database.',flush=True,file=redirect_verbose)
                    go_score=0
                if go_score:
                    score+=go_score
                    score/=2
        if score>1: score=1.0
        if score<0.0001: score=0.0
        return self.get_output(score=score, redirect_verbose=redirect_verbose, console_output=console_output,threshold=threshold,only_return=only_return)





if __name__ == '__main__':
    print(unifunc_folder)
    nlp = UniFunc()
    str1='Responsible for trypanothione reduction'
    str2='Protein associated with trypanothione reduction activity'
    print('Similarity score:',nlp.get_similarity_score(str1,str2,verbose=True,only_return=True))
    str1='Leghemoglobin reductase activity K0002 (EC 0.0.0.0) ID12345 PRK10411.1  '
    str2='Protein associated with trypanothione reductase activity (K0001) ID6789'
    print('Similarity score:',nlp.get_similarity_score(str1,str2,verbose=True,only_return=True))
    print_citation_unifunc()

