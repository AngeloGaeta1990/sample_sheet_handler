import pandas as pd
import Bio
import math
import os
import csv
from ast import While
from cmath import isnan, nan
from configparser import NoOptionError
from msilib import type_binary
from msilib.schema import Class
from pathlib import Path, WindowsPath
from Bio.Seq import Seq




class SampleSheet:
    def __init__(self,input_csv):
        
        # raw dataframe genered from .csv read
        self.raw_dataframe = None
        # daraframe with all rows= NaN filtered
        self.nan_filtered_df = None
        # list of columns name in raw dataframe
        self.columns_list = None
        # list of dictonaris for headers and settings
        self.dictionaries = None 
        # nombers of rows in total dataframe
        self.total_rows = None
        #dictionary containing the header
        self.header_dictionary = None
        # list containint the row index of dataframe in sample sheet
        self.dataframe_rows = [0]
        # position where to find [data] on sample sheet
        self.datasheet_position = None
        # dataframe containing only Data section of sample sheet
        self.data_dataframe = None
        # dataframe v2 
        self.bcl_convert_dataframe_v2 = None
        # list of sample sheet lines
        self.sheet_lines = None




    def read_header_settings(self):
            """
            Created a list of dictonaries for headers and settings
            using samplesheet.csv as input
            """
            #TODO test_input hardcoded replace with -input value
            self.raw_dataframe = pd.read_csv("SampleSheet.csv") 
            self.nan_filtered_df = self.raw_dataframe.dropna(how='all')
            self.total_rows= (len(self.nan_filtered_df))
            self.columns_list = list(self.nan_filtered_df)

    def create_header_dict(self):
                """
                Creates a dictonary for header
                In samplesheet.csv header is in column names rather then rows values-> ad hoc def
                """
                dataframe_dictionaries = []
                header_dict = {}
                if 'Header' in self.columns_list[0]:
                    header_dict['[Header]'] = nan
                    dataframe_start= 0
                    dataframe_end =1
                    
                    usable_indexes = self.dataframe_rows[0:-1]
                    for dataframe in usable_indexes:
                        dataframe_dict= {}
                        for row in range(self.dataframe_rows[dataframe_start], self.dataframe_rows[dataframe_end]):
                            first_column= (self.nan_filtered_df[self.columns_list[0]])
                            first_column_row = first_column.iloc[row]                
                            second_column=((self.nan_filtered_df[self.columns_list[1]]))
                            second_column_row= second_column.iloc[row]
                            dataframe_dict[first_column_row] = second_column_row
                        dataframe_dictionaries.append(dataframe_dict)
                        dataframe_start += 1
                        dataframe_end += 1

                if 'Header' in self.columns_list[0]:
                    first_dict = self.merge_two_dicts(header_dict,dataframe_dictionaries[0])
                    dataframe_dictionaries[0] = first_dict
                
                self.header_dictionary = dataframe_dictionaries
            

    #TODO not really a class function move into miscellanous when ready
    def merge_two_dicts(self,x, y):
        """Given two dictionaries, merge them into a new dict as a shallow copy."""
        z = x.copy()
        z.update(y)
        return z

               

    def find_dataframes(self):
        """
        Finds the row where each dataframe is located
        """
        
        first_column= (self.nan_filtered_df[self.columns_list[0]])
        for row in range(self.total_rows):
            first_column_row = first_column.iloc[row]
            if type(first_column_row) == str:
                # and "Data" not in first_column.iloc[row]
                if "[" in first_column.iloc[row] :
                    self.dataframe_rows.append(row)

    
    def find_samples_data(self):
        """
        convert sample sheet data section into dataframe
        """
        first_column= (self.nan_filtered_df[self.columns_list[0]])
        for row in range(self.total_rows):
            first_column_row = first_column.iloc[row]
            if type(first_column_row) == str:
                if "[" in first_column.iloc[row] and "Data" in first_column.iloc[row] :
                    self.datasheet_position = row
            

    
    def get_sample_data(self):
        """
        convert [data] into dataframe
        """
        datasheet = self.nan_filtered_df.iloc[(self.datasheet_position+1):]
        header_row = datasheet.iloc[[0]]
        data_header = header_row.values
        # creating a 2 element list
        data_header_list= data_header.tolist()
        # only 1 is required
        data_header_list=(data_header_list[0])
        #removing NaN value
        #data_header_list = [item for item in data_header_list if not(pd.isna(item)) == True]
        values_dataframe =(datasheet[1:])
        values_dataframe.columns = data_header_list
        #remove last colums if it is NaN
        if isnan(values_dataframe.columns[-1]) :
           values_dataframe.drop(columns=values_dataframe.columns[-1], axis=1, inplace=True)
        values_dataframe.set_index("Sample_ID",inplace=True)
        self.data_dataframe = values_dataframe
        
        


    def reverse_complement_index2(self):
        """
        takes columns index2 from Data and convert index2 in the reverse complement
        """  
        index2_columns = self.data_dataframe["index2"].tolist()
        reverse_complement_index2 = []
        for index in index2_columns:
            dna_seq = Seq(index)
            rev_comp_index = dna_seq.reverse_complement()
            rev_comp_index = str(rev_comp_index)
            reverse_complement_index2.append(rev_comp_index)
        self.data_dataframe["index2"] = reverse_complement_index2

    def convert_data_to_v2(self):
         """
         conver data section of sample sheet into BCL convert settings
         """
         self.bcl_convert_dataframe_v2 = self.data_dataframe[[ "index","index2"]]


    def turn_headers_into_nested_lists(self):
         """
         created a nested list including the header
         """
         nested_list = [[k, v] for d in self.header_dictionary for k, v in d.items()]
         self.samplesheet_lines = [[x for x in sublist if isinstance(x, str) or (isinstance(x, (int, float)) and not math.isnan(x))] for sublist in nested_list]
         self.samplesheet_lines.append(["[Data]"])
         datasheet_list = [row.tolist() for row in self.data_dataframe.values]
         for row in datasheet_list:
              self.samplesheet_lines.append(row)

    def check_for_illegal_characters(self,sample_name):
        """
        takes a string as input and verifies if there is an illegal characters
        """
        illegal_characters = ["!"," ","#", "$", "%", "&", "'", "(", ")", "*", "+", ",", "-", 
                               ".", "/", ":", ";", ">", "<"]
        for illegal in illegal_characters:
            if illegal in sample_name:
                return True
        return False

    def find_illegal_characters(self):
         """
         browse trough [data] and looks for illegal characters in the sample name
         """
         for index, row in self.data_dataframe.iterrows():
            if self.check_for_illegal_characters(row['Sample_Name']):
                print(row)

            else:
                print("niente")
                
         
        
   


    def rewrite_csv_version2(self, outputfilename):
        """
        this function writes back the .csv files after editing
        """
        self.bcl_convert_dataframe_v2.to_csv(outputfilename)
    
         
    def rewrite_csv_reversed(self,outputfilename):
        """
        this function writes back the .csv files after editing
        """

        #self.data_dataframe.to_csv(outputfilename)
        with open(outputfilename, 'w', newline='') as file:
            writer = csv.writer(file)
            for row in self.samplesheet_lines:
                writer.writerow(row)
 
             
def reverse_complement_index2(filename:"SampleSheet_V1.csv", outputfilename):

    """
    This function converts your index2 in the reverse complement format

    filename = sample sheet filanme type string 
    main function to print sample sheet 
    """
    sample_sheet = SampleSheet(filename) 
    sample_sheet.read_header_settings()
    sample_sheet.find_dataframes()
    sample_sheet.create_header_dict()
    sample_sheet.find_samples_data()
    sample_sheet.get_sample_data()
    sample_sheet.reverse_complement_index2()
    sample_sheet.rewrite_csv_reversed(outputfilename)

def convert_sample_sheet_to_v2(filename:"input V1 sample sheet in .csv format",
                               outputfilename:"filename of your output"):
    """
    This function converts your [Data] section of a V1 SampleSheet into a V2 compatible section
    copy paste this output into your samplesheet

    example commmand
    python main.py convert-sample-sheet-to-v2 my_samplesheet1.csv my_samplesheet2.csv
    my_samplesheet1.csv = input V1 sample sheet in .csv format
    my_samplesheet2.csv = filename of your output

    """
    sample_sheet = SampleSheet(filename) 
    sample_sheet.read_header_settings()
    sample_sheet.find_dataframes()
    sample_sheet.create_header_dict()
    sample_sheet.find_samples_data()
    sample_sheet.get_sample_data()
    sample_sheet.convert_data_to_v2()
    sample_sheet.rewrite_csv_version2(outputfilename)



sample_sheet = SampleSheet("SampleSheet.csv", ) 
sample_sheet.read_header_settings()
sample_sheet.find_dataframes()
sample_sheet.create_header_dict()
sample_sheet.find_samples_data()
sample_sheet.get_sample_data()
sample_sheet.reverse_complement_index2()
sample_sheet.find_illegal_characters()
sample_sheet.turn_headers_into_nested_lists()
sample_sheet.rewrite_csv_reversed("SampleSheet_reversed.csv")






        


