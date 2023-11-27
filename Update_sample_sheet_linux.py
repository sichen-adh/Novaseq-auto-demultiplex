#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import re
import pandas as pd
import glob
import warnings
import sys
#import os
import shutil
import json
import subprocess

def complement_rule(input_sequence):
    complement_rule = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    input_sequence=input_sequence.rstrip()
    complement_sequence=''
    for letter in input_sequence:
        complement_sequence+=complement_rule[letter]
    return complement_sequence

def reverse_complement_rule(input_sequence):
    complement_rule = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'}
    input_sequence=input_sequence.rstrip()
    #print (input_sequence)
    reverse_sequence=input_sequence[::-1]
    reverse_complement_sequence=''

    for letter in reverse_sequence:
        reverse_complement_sequence+=complement_rule[letter]
    return reverse_complement_sequence

def found_common_barcode(i7_index_seq,i5_index_seq,unknown_barcode_df):
    i7_found_df=unknown_barcode_df.loc[unknown_barcode_df['i7_seq']==i7_index_seq]
    i5_found_df=unknown_barcode_df.loc[unknown_barcode_df['i5_seq']==i5_index_seq]
    i7_found_list=i7_found_df.index.values.tolist()
    i5_found_list=i5_found_df.index.values.tolist()
    common_list=list(set(i7_found_list).intersection(i5_found_list))
    #print (common_list)
    return(common_list)


# In[ ]:

############# start to check the demux summary file content #############################
ori_sample_sheet_file=sys.argv[1]
m=re.match('SampleSheet-(.*)-(\d+)nt-(\d+)nt-(\d+)-mismatch.csv',ori_sample_sheet_file)
index1=m.group(2)
index2=m.group(3)
mismatch_info=m.group(4)
run_id=m.group(1)
DemuxSummary_files=glob.glob('./Data/Intensities/BaseCalls/predem_%snt_%snt_%s-mismatch/Stats/DemuxSummary*.txt' %(index1,index2,mismatch_info))
warnings.filterwarnings(action='ignore')

output_tmp_samplesheet_file='./Data/Intensities/BaseCalls/predem_%snt_%snt_%s-mismatch/SampleSheet_tmp.csv' %(index1,index2,mismatch_info)
output_tmp_Samplesheet_f=open(output_tmp_samplesheet_file,'w+')
output_tmp_Samplesheet_f.write('lane,sample_id,i7_index,i5_index,reads_count,comment,\n')

for DemuxSummary_file in DemuxSummary_files:
    match_result=re.match('.*DemuxSummaryF1L(\d+).txt',DemuxSummary_file)
    lane_number=str(match_result.group(1))
    #DemuxSummary_f=open("./%s" %DemuxSummary_file)
    DemuxSummary_f=open(DemuxSummary_file)
    unknown_barcode_df=pd.DataFrame(columns=['i7_seq','i5_seq','reads_count'])
    for line in DemuxSummary_f:
        if line.startswith("### Columns"):
            break

    for line in DemuxSummary_f:
        index_combination,reads_count=line.strip().split('\t')
        #print (index_combination)
        try:
            unknown_barcode_df=unknown_barcode_df.append({'i7_seq':index_combination.split('+')[0],'i5_seq':index_combination.split('+')[1],
                                                     'reads_count':reads_count},ignore_index=True)
        except:
            DemuxSummary_f.close()
            break
        #print (line)

    DemuxSummary_f.close()


    ############################################ collect the barcode info from the sample sheet########################

    ori_sample_barcode_df=pd.DataFrame(columns=['sample_id','i7_index','i5_index'])
    ori_sample_sheet_f=open("./%s" %ori_sample_sheet_file)
    for line in ori_sample_sheet_f:
        if line.startswith('Lane'):
            item=line.split(',')
            #print (item)
            lane_number_index=item.index('Lane')
            sample_id_index=item.index('Sample_ID')
            i7_index_index=item.index('index')
            i5_index_index=item.index('index2')
            #print(sample_id_index,i7_index_index,i5_index_index)
            break

    for line in ori_sample_sheet_f:
        info=line.strip().split(",")
        #print (info[lane_number_index],lane_number)
        if info[lane_number_index]==lane_number:
            ori_sample_barcode_df=ori_sample_barcode_df.append({'sample_id':info[sample_id_index],'i7_index':info[i7_index_index],
                                                            'i5_index':info[i5_index_index]},ignore_index=True)
            #print (line)
    ori_sample_sheet_f.close()



    ########################################## start to match the barcode info ###############################




    for idx,row in ori_sample_barcode_df.iterrows():
        #4 conditions here, original, reverse complement, reverse, complement
        condition=''
        for i7_condition in ['original','reverse complement','reverse','complement']:
            for i5_condition in ['original','reverse complement','reverse','complement']:
                #i7 condition check
                if i7_condition=='original':
                    i7_seq=row['i7_index']
                elif i7_condition=='reverse complement':
                    i7_seq=reverse_complement_rule(row['i7_index'])
                elif i7_condition=='reverse':
                    i7_seq=row['i7_index'][::-1]
                elif i7_condition=='complement':
                    i7_seq=complement_rule(row['i7_index'])

                #i5 condition check
                if i5_condition=='original':
                    i5_seq=row['i5_index']
                elif i5_condition=='reverse complement':
                    i5_seq=reverse_complement_rule(row['i5_index'])
                elif i5_condition=='reverse':
                    i5_seq=row['i5_index'][::-1]
                elif i5_condition=='complement':
                    i5_seq=complement_rule(row['i5_index'])

                if len(found_common_barcode(i7_seq,i5_seq,unknown_barcode_df))!=0:
                    condition='i7 is %s and i5 is %s' %(i7_condition,i5_condition)
                    common_list=found_common_barcode(i7_seq,i5_seq,unknown_barcode_df)

        if condition=='':
            for i7_condition in ['original','reverse complement','reverse','complement']:
                for i5_condition in ['original','reverse complement','reverse','complement']:
                    #i7 condition check
                    if i7_condition=='original':
                        i7_seq=row['i7_index']
                    elif i7_condition=='reverse complement':
                        i7_seq=reverse_complement_rule(row['i7_index'])
                    elif i7_condition=='reverse':
                        i7_seq=row['i7_index'][::-1]
                    elif i7_condition=='complement':
                        i7_seq=complement_rule(row['i7_index'])

                    #i5 condition check
                    if i5_condition=='original':
                        i5_seq=row['i5_index']
                    elif i5_condition=='reverse complement':
                        i5_seq=reverse_complement_rule(row['i5_index'])
                    elif i5_condition=='reverse':
                        i5_seq=row['i5_index'][::-1]
                    elif i5_condition=='complement':
                        i5_seq=complement_rule(row['i5_index'])

                    if len(found_common_barcode(i5_seq,i7_seq,unknown_barcode_df))!=0:
                        condition='i7 is %s and i5 is %s also i5 and i7 are reversed' %(i5_condition,i7_condition)
                        common_list=found_common_barcode(i5_seq,i7_seq,unknown_barcode_df)




            if condition=='':
                print (condition)
                output_tmp_Samplesheet_f.write(','.join([lane_number,row['sample_id'],'','','',condition,'\n']))
            else:
                print (condition)
                output_tmp_Samplesheet_f.write(','.join([lane_number,row['sample_id'],unknown_barcode_df.loc[common_list[0]]['i7_seq'],unknown_barcode_df.loc[common_list[0]]['i5_seq'],
                                              unknown_barcode_df.loc[common_list[0]]['reads_count'],condition,'\n']))
        else:
            print (condition)
            output_tmp_Samplesheet_f.write(','.join([lane_number,row['sample_id'],unknown_barcode_df.loc[common_list[0]]['i7_seq'],unknown_barcode_df.loc[common_list[0]]['i5_seq'],
                                              unknown_barcode_df.loc[common_list[0]]['reads_count'],condition,'\n']))

output_tmp_Samplesheet_f.close()


# In[ ]:


########## Go over the tmp samplesheet file to calculate each comment frequency and the total reads of that comment  #######
output_tmp_Samplesheet_f=open(output_tmp_samplesheet_file)
output_tmp_Samplesheet_f.readline()
#output_sample_sheet_f.readline()
update_sample_sheet_dict={}

for line in output_tmp_Samplesheet_f:
    items=line.split(',')
    lane_number=items[0]
    sample_id=items[1]
    #i7_index=items[2]
    #i5_index=items[3]
    reads_count=items[4]
    comment=items[5].strip()

    sample_id='-'.join(sample_id.split("-")[0:2])
    reads_count=0 if reads_count =='' else int(reads_count)
    comment='original' if comment=='' else comment

    lane_sample_id=lane_number+'+'+sample_id
    if lane_sample_id not in update_sample_sheet_dict.keys():
        update_sample_sheet_dict[lane_sample_id]={}
        update_sample_sheet_dict[lane_sample_id][comment]=[1,reads_count]
    elif comment not in update_sample_sheet_dict[lane_sample_id].keys():
        update_sample_sheet_dict[lane_sample_id][comment]=[1,reads_count]

    else:
        update_sample_sheet_dict[lane_sample_id][comment][0]=update_sample_sheet_dict[lane_sample_id][comment][0]+1
        update_sample_sheet_dict[lane_sample_id][comment][1]=update_sample_sheet_dict[lane_sample_id][comment][1]+reads_count


output_tmp_Samplesheet_f.close()
#print (update_sample_sheet_dict)


# In[ ]:


############ Only keep one comment if there are multiple. Based on the frequency. If there is conexisting of original and other
############ comments, keep the other comment first. It will be double checked in LaneBarcode.html also.
for key,value in update_sample_sheet_dict.items():
    if len(value.values())>1:
        max_frequency=0
        for comment,frequency_reads_list in value.items():
            if comment=='orignal':
                continue
            if frequency_reads_list[0]>max_frequency:
                max_frequency=frequency_reads_list[0]
                kept_comment=comment
                kept_frequency_reads_list=frequency_reads_list
        update_sample_sheet_dict[key]={kept_comment:kept_frequency_reads_list}
#print (update_sample_sheet_dict)


# In[ ]:


#######Xiaoyu script to make the LaneBarcode.html to csv file ###########
def Hcath():
    summary_file=glob.glob("./Data/Intensities/BaseCalls/predem_%snt_%snt_%s-mismatch/Reports/html/*/all/all/all/laneBarcode.html" %(index1,index2,mismatch_info))[0]
    with open("%s" %summary_file, "r") as f:       # read html file
        htmll = f.read()
    html_data = pd.read_html(htmll, match='Sample',header=0, index_col=0)      # match table 2 for lane summary

    for i in html_data:
        table_date = pd.DataFrame(i)
        table_date.to_csv( prj, encoding='utf-8-sig')

def Fmt():
    df = pd.read_csv(prj)
    df.head(2)
    df = df.drop(columns=['Project' , '% of thelane', '% Perfectbarcode', '% One mismatchbarcode', '% PFClusters'], axis=1) # delete columns
    undet = df[df['Sample'] == 'Undetermined'].index.tolist()
    df = df.drop(undet) # delete undetermined rows
	#multiple lanes sum up
    #aggregation_functions = {'Lane':'first','Barcode sequence': 'first', 'PF Clusters': 'sum', 'Yield (Mbases)': 'sum', '% >= Q30bases': 'first', 'Mean QualityScore': 'first'}
    #df_new = df.groupby('Sample',as_index=False).aggregate(aggregation_functions).reindex(columns=df.columns)
	#write to csv
    df.to_csv(prj, encoding='utf-8-sig', index=False)

prj='./Data/Intensities/BaseCalls/predem_%snt_%snt_%s-mismatch/QCsummary.csv' %(index1,index2,mismatch_info)
#load html report
Hcath()
#format csv
Fmt()





######## To compare the reads number of the updated samples to the orignial ones. If the updated ones less than the original ones, no need to update samplesheet
# first to get the reads number that groupy by lane and project_id
QC_summary_f=open('./Data/Intensities/BaseCalls/predem_%snt_%snt_%s-mismatch/QCsummary.csv' %(index1,index2,mismatch_info))
QC_summary_f.readline()
reads_count_dict={}
for line in QC_summary_f:
    items=line.split(',')
    lane_num=items[0]
    sample_id=items[1]
    sample_id='-'.join(sample_id.split("-")[0:2])
    reads_count=int(items[3])

    lane_sample_id=lane_num+'+'+sample_id
    if lane_sample_id not in reads_count_dict.keys():
        reads_count_dict[lane_sample_id]=reads_count
    elif lane_sample_id in reads_count_dict.keys():
        reads_count_dict[lane_sample_id]+=reads_count

QC_summary_f.close()
#if the orignal barcode got more reads, then change the comment back to original.
for lane_sample_id,updated_value in update_sample_sheet_dict.items():
    #print (reads_count_dict[lane_sample_id],list(updated_value.values()))
    if reads_count_dict[lane_sample_id]>list(updated_value.values())[0][1]:
        update_sample_sheet_dict[lane_sample_id]={'original':list(updated_value.values())[0]}


print (update_sample_sheet_dict)
if len(update_sample_sheet_dict)==0:
    quit()
with open('./Data/Intensities/BaseCalls/predem_%snt_%snt_%s-mismatch/Update_sample_sheet_dict.txt' %(index1,index2,mismatch_info),'w+') as update_sample_sheet_dict_f:
    for key,value in update_sample_sheet_dict.items():
        update_sample_sheet_dict_f.write(key)
        update_sample_sheet_dict_f.write('\t')
        update_sample_sheet_dict_f.write(json.dumps(value))
        update_sample_sheet_dict_f.write('\n')
        print (json.dumps(value))




########### modify the samplesheet ##############
ori_sample_sheet_f=open("./%s" %ori_sample_sheet_file)
update_sample_sheet_file="./Data/Intensities/BaseCalls/predem_%snt_%snt_%s-mismatch/Update_%s" %(index1,index2,mismatch_info,ori_sample_sheet_file)
update_sample_sheet_f=open(update_sample_sheet_file,'w+')

for line in ori_sample_sheet_f:
    if not line.startswith('Lane'):
        update_sample_sheet_f.write(line)
        continue
    elif line.startswith('Lane'):
        items=line.split(',')
        Lane_index=items.index('Lane')
        Sample_ID_index=items.index('Sample_ID')
        i7_seq_index=items.index('index')
        i5_seq_index=items.index('index2')
        update_sample_sheet_f.write(line)
        break

for line in ori_sample_sheet_f:
    items=line.split(',')
    sample_id=items[Sample_ID_index]
    sample_id='-'.join(sample_id.split("-")[0:2])
    lane_sample_id=items[Lane_index]+'+'+sample_id
    #print (list[update_sample_sheet_dict[lane_sample_id].keys()][0][0])
    comment_list=[key for key in update_sample_sheet_dict[lane_sample_id].keys()]
    comment=comment_list[0]
    print (comment)

    if comment=='original':
        update_sample_sheet_f.write(line)

    else:
        condition1=comment.split(',')[0]
        m=re.match('i7 is (.*) and i5 is (.*)',condition1)
        i7_condition=m.group(1)
        i5_condition=m.group(2)

        if i7_condition=='original':
            i7_update_seq=items[i7_seq_index]
        elif i7_condition=='reverse':
            i7_update_seq=items[i7_seq_index][::-1]
        elif i7_condition=='reverse complement':
            i7_update_seq=reverse_complement_rule(items[i7_seq_index])
        elif i7_condition=='complement':
            i7_update_seq=complement_rule(items[i7_seq_index])

        if i5_condition=='original':
            i5_update_seq=items[i5_seq_index]
        elif i5_condition=='reverse':
            i5_update_seq=items[i5_seq_index][::-1]
        elif i5_condition=='reverse complement':
            i5_update_seq=reverse_complement_rule(items[i5_seq_index])
        elif i7_condition=='complement':
            i5_update_seq=complement_rule(items[i5_seq_index])

        try:
            condition2=comment.split(',')[1]
            items[i7_seq_index]=i5_update_seq
            items[i5_seq_index]=i7_update_seq
        except:
            items[i7_seq_index]=i7_update_seq
            items[i5_seq_index]=i5_update_seq
        update_sample_sheet_f.write(','.join(items))

ori_sample_sheet_f.close()
update_sample_sheet_f.close()

#delete the original samplesheet, mv the updated sample sheet here.
shutil.move(update_sample_sheet_file,ori_sample_sheet_file)
send_report='mailx -A ./Data/Intensities/BaseCalls/predem_%snt_%snt_%s-mismatch/Update_sample_sheet_dict.txt -A ./Data/Intensities/BaseCalls/predem_%snt_%snt_%s-mismatch/SampleSheet_tmp.csv -s %s_predem_%snt_%snt_%s-mismatch-update si.chen@admerahealth.com shuang.wu@admerahealth.com li.mengxiang@admerahealth.com </dev/null' %(index1,index2,mismatch_info,index1,index2,mismatch_info,run_id,index1,index2,mismatch_info)
subprocess.check_call(send_report,shell=True)
