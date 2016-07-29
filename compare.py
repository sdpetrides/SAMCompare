# SAMCompare
# @author: Stephen Petrides

import sys

class SAMRecord:
	"""Holds necessary information from SAM record"""
	def __init__(self, FLAG, RNAME, POS, MAPQ, SAM_file, AS='', NM=''):
		self.FLAG = FLAG
		self.MAPQ = MAPQ
		self.RNAME = RNAME
		if len(POS) < 8:
			self.POS = POS + '\t'
		else:
			self.POS = POS
		if AS != '':
			self.AS = 'AS' + AS[4:]
			if SAM_file == 1:
				stats_dict_1['Alignments']+=1
			else:
				stats_dict_2['Alignments']+=1
		else:
			self.AS = AS
		if NM != '':	
			for field in NM:
				if field.startswith('NM'):
					self.NM = field
					break
		else:
			self.NM = NM
		if SAM_file == 1:
			stats_dict_1['Records']+=1
		else:
			stats_dict_2['Records']+=1

	def __str__(self):
		return(
				self.FLAG + '\t' + 
				self.MAPQ + '\t' +
				self.RNAME + '\t' + 
				self.POS + '\t' +
				self.AS + '\t' + 
				self.NM
			)

	def __eq__(self, other):
		if self.RNAME != other.RNAME:
			return False
		elif self.POS != other.POS:
			return False
		else:
			return True

class Header:
	"""Holds paired lists with records from a SAM file"""

	def __init__(self, header_name):
		self.header_name = header_name
		self.first_ip_list = []
		self.second_ip_list = []

	def add_record(self, SAMRecord, in_pair):
		if in_pair == 1:
			self.first_ip_list.append(SAMRecord)
		else:
			self.second_ip_list.append(SAMRecord)

	def __str__(self):
		string = ''
		while self.first_ip_list or self.second_ip_list:
			try:
				string+= '   ' + str(self.first_ip_list.pop(0)) + '\n'
			except:
				pass
			try:
				string+= '   ' + str(self.second_ip_list.pop(0)) + '\n'
			except:
				pass
		return self.header_name + '\n' + string

	def __eq__(self, other):
		if self.header_name != other.header_name:
			return False

		elif (len(self.first_ip_list) != len(other.first_ip_list) or 
				len(self.second_ip_list) != len(other.second_ip_list)):
			diff_dict['Number of Alignments']+=1
			return False
		
		match_count = 0
		match_threshold = len(self.first_ip_list)
		for record_self in self.first_ip_list:
			for record_other in other.first_ip_list:
				if record_self == record_other:
					match_count+=1

		if match_count < match_threshold:
			diff_dict['Mapped Position']+=1
			return False

		match_count = 0
		match_threshold = len(self.second_ip_list)	
		for record_self in self.second_ip_list:
			for record_other in other.second_ip_list:
				if record_self == record_other:
					match_count+=1

		if match_count < match_threshold:
			diff_dict['Mapped Position']+=1
			return False

		diff_dict['No Change']+=1

		return True

def SAM_rec_to_SAM_obj(header, fields, FLAG_int, SAM_file):
	try:
		if FLAG_int & 64:
			header.add_record(
				SAMecord(
					str(bin(FLAG_int)), 
					fields[2], 
					fields[3], 
					fields[4], 
					SAM_file, 
					AS=fields[11], 
					NM=(fields[16], 
					fields[17],
					),
				), 
				1,
			)
		else:
			header.add_record(
				SAMRecord(
					str(bin(FLAG_int)),
					fields[2], 
					fields[3], 
					fields[4], 
					SAM_file, 
					AS=fields[11], 
					NM=(fields[16], 
						fields[17],
					),
				), 
				2,
			)
	except:
		if FLAG_int & 64:
			header.add_record(
				SAMRecord(
					str(bin(FLAG_int)), 
					fields[2], 
					fields[3], 
					fields[4], 
					SAM_file,
				), 
				1,
			)
		else:
			header.add_record(
				SAMRecord(
					str(bin(FLAG_int)),
					fields[2], 
					fields[3], 
					fields[4], 
					SAM_file,
				), 
				2,
			)

def calc_print_stats():
	names_1 = [
		'Headers',
		'Records',
		'Records per Header',
		'Alignments',
		'Alignments per Header',
		'Unmapped',
	]

	names_2 = [
		'No Change',
		'Number of Alignments',
		'Mapped Position',
	]

	stats_dict_1["Alignments per Header"] = str(
		format(
			stats_dict_1["Alignments"] / stats_dict_1["Headers"],
			'.3f',
			)
		)
	stats_dict_2["Alignments per Header"] = str(
		format(
			stats_dict_2["Alignments"] / stats_dict_2["Headers"],
			'.3f',
			)
		)
	stats_dict_1["Records per Header"] = str(
		format(
			stats_dict_1["Records"] / stats_dict_1["Headers"],
			'.3f',
			)
		)
	stats_dict_2["Records per Header"] = str(
		format(
			stats_dict_2["Records"] / stats_dict_2["Headers"],
			'.3f',
			) 
		)
	stats_dict_1["Unmapped"] = str(
		stats_dict_1["Records"] - stats_dict_1["Alignments"]
		)
	stats_dict_2["Unmapped"] = str(
		stats_dict_2["Records"] - stats_dict_2["Alignments"]
		)

	print('\nSAM FILE\t\t\t1\t2')
	for key in names_1:
		if len(key) < 8:
			print(
				key + '\t\t\t\t' + 
				str(stats_dict_1[key]) + '\t' + 
				str(stats_dict_2[key])
			)
		elif len(key) < 16:
			print(
				key + '\t\t\t' + 
				str(stats_dict_1[key]) + '\t' + 
				str(stats_dict_2[key])
			)
		else: 
			print(
				key + '\t\t' + 
				str(stats_dict_1[key]) + '\t' + 
				str(stats_dict_2[key])
			)
	print('\nDIFFERENCE\t\t#\t%')
	for key in names_2:
		if len(key) < 8:
			print(
					key 
					+ '\t\t\t' 
					+ str(diff_dict[key]) 
					+ '\t' 
					+ str(
						format(
							100*(diff_dict[key]/stats_dict_1['Headers']),
						'.2f',
						)
					) 
					+ '\t' 
					+ str(
						format(
							100*(diff_dict[key]/stats_dict_2['Headers']),
						'.2f',
						)
					)
				)
		elif len(key) < 16:
			print(
					key 
					+ '\t\t' 
					+ str(diff_dict[key]) 
					+ '\t' 
					+ str(
						format(
							100 * (diff_dict[key] / stats_dict_1['Headers']), 
							'.2f',
						)
					)
					+ '\t' 
					+ str(
						format(
							100 * (diff_dict[key] / stats_dict_2['Headers']),
							'.2f',
						)
					)
				)
		else:
			print(
					key 
					+ '\t' 
					+ str(diff_dict[key]) 
					+ '\t' 
					+ str(
						format(
							100 * (diff_dict[key] / stats_dict_1['Headers']),
							'.2f',
						)
					) 
					+ '\t' 
					+ str(
						format(
							100 * (diff_dict[key] / stats_dict_2['Headers']),
							'.2f',
						)
					)
				)

def open_SAM_file(SAM_str):
	if SAM_str[-4:] != '.sam':
		print(	
			"Not valid SAM file. Must have .sam extension.\n" +
			"Error: "
			+ SAM_str
		)
		sys.exit(0)
	try:
		SAM_file = open(SAM_str, 'rt')
	except:
		print(	
			"Could not open SAM file.\n" +
			"Error: "
			+ SAM_str
		)
		sys.exit(0)
	return SAM_file


stats_dict_1 = {
	'Headers' : 0,
	'Alignments' : 0,
	'Alignments per Header' : '',
	'Records per Header' : '',
	'Records' : 0,
	'Unmapped' : 0,
}

stats_dict_2 = {
	'Headers' : 0,
	'Alignments' : 0,
	'Alignments per Header' : '',
	'Records per Header' : '',
	'Records' : 0,
	'Unmapped' : 0,
}

diff_dict = {
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

	manual = [
		"\nCommands and Options", 
		"\n\n    Basic Command Structure:", 
		"\n\n\t$ python3 compare.py <flag> <path>",
		"\n\t\t# Flags can come in any order but must be followed by ", 
		"corresponding information",
		"\n\t\t# Running the program without any flags or ", 
		"additional arguments will print the manual",
		"\n    Flags", 
		"\n\t-man", 
		"\n\t\t# Prints Commands and Options section of MANUAL.", 
		"\n\t-1 <path> (required)",
		"\n\t\t# Followed by the path/to/file_name_1.sam", 
		"\n\t-2 <path> (required)",
		"\n\t\t# Followed by the path/to/file_name_2.sam", 
		"\n\t-ex_count <int> (optional)"
		"\n\t\t# Followed by the number of example alignments to print",
		"\n\t\t# Default is 0 example alignments to print", 
		"\n\t-no_stats (optional)",
		"\n\t\t# Does not print statistics"
	]

	while param_count <= argc-1:
		param_str = sys.argv[param_count]
		param_set.add(param_str)
		param_count +=1
		if param_str[0] is not '-':
			print(
				"Parameters not entered correctly. " +
				"See MANUAL for proper syntax.\n" +
				"Error: "
				 + param_str
			)
			sys.exit(0)
		else:
			if param_str == '-man':
				print(''.join(manual))
				sys.exit(0)
			elif param_str == '-1':
				param_str = sys.argv[param_count]
				param_count +=1
				SAM_1_str = param_str
			elif param_str == '-2':
				param_str = sys.argv[param_count]
				param_count +=1
				SAM_2_str = param_str
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
				print(
					"Parameters not entered correctly. " +
					"See MANUAL for proper syntax.\n" +
					"Error: "
					+ param_str
				)
				sys.exit(0)

	if not param_set:
		print(''.join(manual))
		sys.exit(0)

	master_dict_1 = {}

	SAM_file_1 = open_SAM_file(SAM_1_str)

	for line in SAM_file_1:
		fields = line.split()
		first_word = fields[0]
		if first_word[0] != '@':
			FLAG_int = int(fields[1])

			header = Header(fields[0])
			if header.header_name in master_dict_1:
				header = master_dict_1[header.header_name]
				SAM_rec_to_SAM_obj(
					header, 
					fields, 
					FLAG_int,
					1,
					)
			else:
				SAM_rec_to_SAM_obj(
					header, 
					fields, 
					FLAG_int,
					1,
					)
				master_dict_1[header.header_name] = header
				stats_dict_1['Headers']+=1
	SAM_file_1.close()

	SAM_file_2 = open_SAM_file(SAM_2_str)
	header_count = 0
	header = Header('')
	prev_header = Header('')

	for line in SAM_file_2:
		fields = line.split()
		first_word = fields[0]
		if first_word[0] != '@':
			FLAG_int = int(fields[1])

			if (prev_header.header_name == '' 
				or prev_header.header_name != fields[0]):
				header = Header(fields[0])
				stats_dict_2['Headers']+=1

				SAM_rec_to_SAM_obj(
					header, 
					fields, 
					FLAG_int,
					2,
					)

				if prev_header.header_name != '':
					header_1 = master_dict_1.pop(prev_header.header_name)
					if header_1 == prev_header:
						pass
					elif header_count < header_limit:
						print(header_1)
						print(prev_header)
						header_count+=1	
			else:
				SAM_rec_to_SAM_obj(
					header, 
					fields, 
					FLAG_int,
					2,
					)
			prev_header = header

	header_1 = master_dict_1.pop(prev_header.header_name)
	if header_1 == prev_header:
		pass
	elif header_count < header_limit:
		print(header_1)
		print(prev_header)
		header_count+=1	
	
	SAM_file_2.close()

	if not no_stats:
		calc_print_stats()

if __name__ == "__main__":
	main()