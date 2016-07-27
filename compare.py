# SAMCompare
# @author: Stephen Petrides

import sys
import random

class SAM_RECORD_OBJ:
	"""Holds important information from SAM record"""
	def __init__(self, FLAG, RNAME, POS, MAPQ, sam_file, AS='', NM='', ):
		self.FLAG = FLAG
		self.MAPQ = MAPQ
		self.RNAME = RNAME
		if len(POS) < 8:
			self.POS = POS + '\t'
		else:
			self.POS = POS
		if AS != '':
			self.AS = 'AS' + AS[4:]
			if sam_file == 1:
				STATS_1['Alignments']+=1
			else:
				STATS_2['Alignments']+=1
		else:
			self.AS = AS
		if NM != '':	
			for field in NM:
				if field.startswith('NM'):
					self.NM = field
					break
		else:
			self.NM = NM
		if sam_file == 1:
			STATS_1['Records']+=1
		else:
			STATS_2['Records']+=1

	def __str__(self):
		return self.FLAG + '\t' + self.MAPQ + '\t' + self.RNAME + '\t' + self.POS + '\t' + self.AS + '\t' + self.NM

	def __eq__(self, other):
		if self.RNAME != other.RNAME:
			return False
		elif self.POS != other.POS:
			return False
		else:
			return True

class HEADER_OBJ:
	"""Holds two lists with records from a SAM file"""

	def __init__(self, header_name):
		self.header_name = header_name
		self.first_ip_list = []
		self.second_ip_list = []

	def add_record(self, SAM_RECORD_OBJ, in_pair):
		if in_pair == 1:
			self.first_ip_list.append(SAM_RECORD_OBJ)
		else:
			self.second_ip_list.append(SAM_RECORD_OBJ)

	def __str__(self):
		string = ''
		while self.first_ip_list or self.second_ip_list:
			try:
				string = string + '   ' + str(self.first_ip_list.pop(0)) + '\n'
			except:
				pass
			try:
				string = string + '   ' + str(self.second_ip_list.pop(0)) + '\n'
			except:
				pass
		return self.header_name + '\n' + string

	def __eq__(self, other):
		if self.header_name != other.header_name:
			return False

		elif len(self.first_ip_list) != len(other.first_ip_list) or len(self.second_ip_list) != len(other.second_ip_list):
			DIFF_STATS['Number of Alignments']+=1
			return False
		
		match_count = 0
		match_threshold = len(self.first_ip_list)
		for record_self in self.first_ip_list:
			for record_other in other.first_ip_list:
				if record_self == record_other:
					match_count+=1

		if match_count < match_threshold:
			DIFF_STATS['Mapped Position']+=1
			return False

		match_count = 0
		match_threshold = len(self.second_ip_list)	
		for record_self in self.second_ip_list:
			for record_other in other.second_ip_list:
				if record_self == record_other:
					match_count+=1

		if match_count < match_threshold:
			DIFF_STATS['Mapped Position']+=1
			return False

		DIFF_STATS['No Change']+=1
		return True

def calc_print_STATS():
	STATS_1["Alignments per Header"] = str(STATS_1["Alignments"]/STATS_1["Headers"])
	STATS_2["Alignments per Header"] = str(STATS_2["Alignments"]/STATS_2["Headers"])
	STATS_1["Records per Header"] = str(STATS_1["Records"]/STATS_1["Headers"])
	STATS_2["Records per Header"] = str(STATS_2["Records"]/STATS_2["Headers"])
	STATS_1["Unmapped"] = str(STATS_1["Records"]-STATS_1["Alignments"])
	STATS_2["Unmapped"] = str(STATS_2["Records"]-STATS_2["Alignments"])

	print('\nSAM FILE\t\t\t1\t2')
	for key in STATS_1:
		if len(key) < 8:
			print(key + '\t\t\t\t' + str(STATS_1[key]) + '\t' + str(STATS_2[key]))
		elif len(key) < 16:
			print(key + '\t\t\t' + str(STATS_1[key]) + '\t' + str(STATS_2[key]))
		else: 
			print(key + '\t\t' + str(STATS_1[key]) + '\t' + str(STATS_2[key]))
	print('\nDIFFERENCE\t\t#\t%')
	for key in DIFF_STATS:
		if len(key) < 8:
			print(key + '\t\t\t' + str(DIFF_STATS[key]) + '\t' + str(100*(DIFF_STATS[key]/STATS_1['Headers'])) + '\t' + str(100*(DIFF_STATS[key]/STATS_2['Headers'])))
		elif len(key) < 16:
			print(key + '\t\t' + str(DIFF_STATS[key]) + '\t' + str(100*(DIFF_STATS[key]/STATS_1['Headers'])) + '\t' + str(100*(DIFF_STATS[key]/STATS_2['Headers'])))
		else:
			print(key + '\t' + str(DIFF_STATS[key]) + '\t' + str(100*(DIFF_STATS[key]/STATS_1['Headers'])) + '\t' + str(100*(DIFF_STATS[key]/STATS_2['Headers'])))

STATS_1 = {
	'Headers' : 0,
	'Alignments' : 0,
	'Alignments per Header' : '',
	'Records per Header' : '',
	'Records' : 0,
	'Unmapped' : 0,
}

STATS_2 = {
	'Headers' : 0,
	'Alignments' : 0,
	'Alignments per Header' : '',
	'Records per Header' : '',
	'Records' : 0,
	'Unmapped' : 0,
}

DIFF_STATS = {
	'Number of Alignments' : 0,
	'Mapped Position' : 0,
	'No Change' : 0,
}

def main():

	argc = len(sys.argv)
	param_set = set()
	param_count = 1
	header_limit = 0
	no_stats = 0

	manual = ["\nCommands and Options", "\n\n\tBasic Command Structure:", "\n\n\t\t$ python3 compare.py <flag> <path>",
	"\n\t\t\t# Flags can come in any order but must be followed by corresponding information",
	"\n\tFlags", "\n\t\t-man", "\n\t\t\t# Prints Commands and Options section of MANUAL.", "\n\t\t-1 <path> (required)",
	"\n\t\t\t# Followed by the path/to/file_name_1.sam", "\n\t\t-2 <path> (required)",
	"\n\t\t\t# Followed by the path/to/file_name_2.sam", "\n\t\t-ex_count <int> (optional)"
	"\n\t\t\t# Followed by the number of example alignments to print",
	"\n\t\t\t# Default is 0 example alignments to print", "\n\t\t-no_stats (optional)",
	"\n\t\t\t# Does not print statistics"]

	while param_count <= argc-1:
		param_str = sys.argv[param_count]
		param_set.add(param_str)
		param_count +=1
		if param_str[0] is not '-':
			print("Parameters not entered correctly. See MANUAL for proper syntax.\nError: " + param_str)
			sys.exit(0)
		else:
			if param_str == '-man':
				print(''.join(manual))
				sys.exit(0)
			elif param_str == '-1':
				param_str = sys.argv[param_count]
				param_count +=1
				sam_1_str = param_str
			elif param_str == '-2':
				param_str = sys.argv[param_count]
				param_count +=1
				sam_2_str = param_str
			elif param_str == '-header':
				param_str = sys.argv[param_count]
				param_count +=1
				header_str = param_str
			elif param_str == '-ex_count':
				param_str = sys.argv[param_count]
				param_count +=1
				header_limit = int(param_str)
			elif param_str == '-no_stats':
				no_stats = 1
			else:
				print("Parameters not entered correctly. See MANUAL for proper syntax.\nError: " + param_str)
				sys.exit(0)

	MASTER_DICT_1 = {}
	MASTER_DICT_2 = {}

	sam_file_1 = open(sam_1_str, 'rt')
	for line in sam_file_1:
		fields = line.split()
		first_word = fields[0]
		if first_word[0] != '@':
			header = HEADER_OBJ(fields[0])
			if header.header_name in MASTER_DICT_1:
				header = MASTER_DICT_1[header.header_name]
				FLAG_bin = str(bin(int(fields[1])))
				try:
					if FLAG_bin[-7] == '1':
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 1, AS=fields[11], NM=(fields[16], fields[17])), 1)
					else:
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 1, AS=fields[11], NM=(fields[16], fields[17])), 2)
				except:
					if FLAG_bin[-7] == '1':
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 1), 1)
					else:
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 1), 2)
				MASTER_DICT_1[header.header_name] = header
			else:
				FLAG_bin = str(bin(int(fields[1])))
				try:
					if FLAG_bin[-7] == '1':
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 1, AS=fields[11], NM=(fields[16], fields[17])), 1)
					else:
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 1, AS=fields[11], NM=(fields[16], fields[17])), 2)
				except:
					if FLAG_bin[-7] == '1':
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 1), 1)
					else:
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 1), 2)
				MASTER_DICT_1[header.header_name] = header
				STATS_1['Headers']+=1
	sam_file_1.close()

	sam_file_2 = open(sam_2_str, 'rt')

	for line in sam_file_2:
		fields = line.split()
		first_word = fields[0]
		if first_word[0] != '@':
			header = HEADER_OBJ(fields[0])
			if header.header_name in MASTER_DICT_2:
				header = MASTER_DICT_2[header.header_name]
				FLAG_bin = str(bin(int(fields[1])))
				try:
					if FLAG_bin[-7] == '1':
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 2, AS=fields[11], NM=(fields[16], fields[17])), 1)
					else:
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 2, AS=fields[11], NM=(fields[16], fields[17])), 2)
				except:
					if FLAG_bin[-7] == '1':
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 2), 1)
					else:
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 2), 2)
				MASTER_DICT_2[header.header_name] = header
			else:
				FLAG_bin = str(bin(int(fields[1])))
				try:
					if FLAG_bin[-7] == '1':
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 2, AS=fields[11], NM=(fields[16], fields[17])), 1)
					else:
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 2, AS=fields[11], NM=(fields[16], fields[17])), 2)
				except:
					if FLAG_bin[-7] == '1':
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 2), 1)
					else:
						header.add_record(SAM_RECORD_OBJ(FLAG_bin, fields[2], fields[3], fields[4], 2), 2)
				MASTER_DICT_2[header.header_name] = header
				STATS_2['Headers']+=1
	sam_file_2.close()

	header_count = 0
	keys_list = list(MASTER_DICT_1.keys())
	while MASTER_DICT_1:
		random_key = random.choice(keys_list)
		keys_list.remove(random_key)
		entry_1 = MASTER_DICT_1.pop(random_key)
		entry_2 = MASTER_DICT_2.pop(random_key)
		if entry_1 == entry_2:
			continue
		elif header_count < header_limit:
			print(entry_1)
			print(entry_2)
			header_count+=1

	if not no_stats:
		calc_print_STATS()

if __name__ == "__main__":
	main()