#!/usr/bin/python3
# coding: utf-8

# In[1]:


def nova_predem():
    subprocess.check_call('ls',shell=True)
    global error_happened
    error_happened=None
    tree=ET.parse('RunInfo.xml')
    root=tree.getroot()
    for child in root:
        print(child.tag,child.attrib)

    #to check if it read the index so that we can do demulitplex
    length_can_predem=0
    for reads in root.iter('Read'):
        reads_dict=reads.attrib
        Number=int(reads_dict['Number'])
        NumCycles=int(reads_dict['NumCycles'])
        IsIndexedRead=reads_dict['IsIndexedRead']
        if Number==1 or IsIndexedRead=='Y':
            length_can_predem+=NumCycles
    print ("start to do predem after %d cycles" %length_can_predem)

    bcl_folder='./Data/Intensities/BaseCalls/L001/'

    current_pos=0
    for item in os.scandir(bcl_folder):
        if item.is_dir():
            #print(item.path)
            m=re.match('.*/C(\d+).1',item.path)
            if int(m.group(1))> current_pos:
                current_pos =int(m.group(1))
    print ("current %d cycles" %current_pos)

    if current_pos> length_can_predem:
        bcl_copy_complete=check_missing_bcl(length_can_predem)
        if bcl_copy_complete==False:
            return
        error_happened=False
        print ('we can do predem!')
    # check done

        sample_sheets=glob.glob("SampleSheet*mismatch.csv")
        current_path=os.getcwd()
        print (current_path)
        try:
            flowcell_ID=re.split(r'[-_]',current_path)[-1][1:]
        except:
            print ('samplesheet is wrong')
            error_happened=True
            return

        for item in sample_sheets:
            #the error_happened is defined false before the demutiplex start
            print (item)
            m=re.match('SampleSheet-(.*)-(\d+)nt-(\d+)nt-(\d+)-mismatch.csv',item)
            print (m.group(2),m.group(3))
            run_name=m.group(1)
            index1=m.group(2)
            index2=m.group(3)
            mismatch_info=m.group(4)

#################################### change the bcl2fastq parameter here #####################
            if index1=='19' and index2=='10':
                bases_mask="y1n*,i10n*,i10n*,n*"

            elif index1=='5' and index2=='8':
                bases_mask="i5y1n*,i8n*,n*,n*"

            elif index1=='0':
                bases_mask="y*,n*,i%sn*,n*" %index2

            elif index2 !='0':
                bases_mask="y1n*,i%sn*,i%sn*,n*" %(index1,index2)

            elif index2 =='0':
                bases_mask="y1n*,i%sn*,n*,n*" %(index1)
################################# finshed change the parameter ###################################
            try:
                bcl2fastq='bcl2fastq --sample-sheet %s --use-bases-mask %s --barcode-mismatch %s -o Data/Intensities/BaseCalls/predem_%snt_%snt_%s-mismatch 2>>%s_%s_%s-mismatch-predem.log' %(item,bases_mask,mismatch_info,index1,index2,mismatch_info,index1,index2,mismatch_info)
                print (bcl2fastq)
                subprocess.check_call(bcl2fastq,shell=True)
            except:
                error_happened=True
                continue
            send_report='mailx -A Data/Intensities/BaseCalls/predem_{}nt_{}nt_{}-mismatch/Reports/html/*{}/all/all/all/laneBarcode.html $(printf -- "-A %s " Data/Intensities/BaseCalls/predem_{}nt_{}nt_{}-mismatch/Stats/DemuxSummary*.txt) -s {}_predem_{}nt_{}nt_{}-mismatch si.chen@admerahealth.com  shuang.wu@admerahealth.com li.mengxiang@admerahealth.com</dev/null'.format(index1,index2,mismatch_info,flowcell_ID,index1,index2,mismatch_info,run_name,index1,index2,mismatch_info)
            subprocess.check_call(send_report,shell=True)
            print (send_report)

            update_samplesheet='python3 /mnt/novaoutput2/si_command/Update_sample_sheet_linux.py %s' %item
            subprocess.check_call(update_samplesheet,shell=True)


        def delete_predem_fastq(pattern,path):
            for root,dirs,files in os.walk(path):
                for name in files:
                    if fnmatch.fnmatch(name,pattern):
                        file_path=os.path.join(root,name)
                        print ("deleting %s" %file_path)
                        os.remove(file_path)

        file_name_pattern='*fastq.gz'
        current_path=os.getcwd()
        start_path='%s/Data/Intensities/BaseCalls/' %current_path
        delete_predem_fastq(file_name_pattern,start_path)



def nova_dem():
    global error_happened
    error_happened=None
    tree=ET.parse('RunInfo.xml')
    root=tree.getroot()
    for child in root:
        print(child.tag,child.attrib)

    #to check if it read the index so that we can do demulitplex
    length_can_dem=0
    for reads in root.iter('Read'):
        reads_dict=reads.attrib
        Number=int(reads_dict['Number'])
        NumCycles=int(reads_dict['NumCycles'])
        IsIndexedRead=reads_dict['IsIndexedRead']
        length_can_dem+=NumCycles
    print ('total %d cycles' %length_can_dem)

    bcl_folder='./Data/Intensities/BaseCalls/L001/'

    current_pos=0
    for item in os.scandir(bcl_folder):
        if item.is_dir():
            #print(item.path)
            m=re.match('.*/C(\d+).1',item.path)
            if int(m.group(1))> current_pos:
                current_pos =int(m.group(1))
    print ("Current %d cycle" %current_pos)

    if current_pos== length_can_dem:
        bcl_copy_complete=check_missing_bcl(length_can_dem)
        if bcl_copy_complete==False:
            return
        print ('we can do dem!')
        #The rest is to run bcl2fastq
        sample_sheets=glob.glob("SampleSheet*mismatch.csv")
        current_path=os.getcwd()
        print (current_path)
        try:
            flowcell_ID=re.split(r'[-_]',current_path)[-1][1:]
        except:
            error_happened=True
            print ('sample sheet is wrong')
            return
        ##if bcl files not completed, then need to wait to demutiplex
        sample_sheet_first=sample_sheets[0]
        #if 'XP' in sample_sheet_first:
        #    if not os.path.exists('CopyComplete.txt'):
        #        return

        error_happened=False
        for item in sample_sheets:
            #print (item)
            sample_sheet_name=item.split(".")[0]
            m=re.match('SampleSheet-(.*)-(\d+)nt-(\d+)nt-(\d+)-mismatch.csv',item)
            print (m.group(2),m.group(3))
            run_name=m.group(1)
            index1=m.group(2)
            index2=m.group(3)
            mismatch_info=m.group(4)

            ########check the parameter here############
            if  index1=='19' and index2=='10':
                bases_mask='y*,i10y9n*,i10n*,y*'

            elif  index1=='5' and index2=='8':
                bases_mask='i5y*,i8n*,n*,y*'

            elif index1 =='0':
                bases_mask='y*,n*,i%sn*,y*' %index2

            elif index2 !='0':
                bases_mask='y*,i%sn*,i%sn*,y*' %(index1,index2)

            elif index2=='0':
                bases_mask='y*,i%sn*,n*,y*' %index1
        ##################################################
            try:
                bcl2fastq='bcl2fastq --sample-sheet %s --use-bases-mask %s -p 40 --barcode-mismatches %s -o /mnt/novaoutput2/DeliverTeam_Output/6-Dmx-Fastq/%s 2>>%s_%s_%s-mismatch_dem_.log' %(item,bases_mask,mismatch_info,sample_sheet_name,index1,index2,mismatch_info)
                print (bcl2fastq)
                subprocess.check_call(bcl2fastq,shell=True)

            except:
                error_happened=True
                continue

            send_report='mailx -A /mnt/novaoutput2/DeliverTeam_Output/6-Dmx-Fastq/{}/Reports/html/*{}/all/all/all/laneBarcode.html $(printf -- "-A %s " /mnt/novaoutput2/DeliverTeam_Output/6-Dmx-Fastq/{}/Stats/DemuxSummary*.txt) -s {}_demultiplex_result si.chen@admerahealth.com shuang.wu@admerahealth.com jingling.hou@admerahealth.com yaping.feng@admerahealth.com kevin.sun@admerahealth.com haixin.shu@admerahealth.com yaoqi.li@admerahealth.com li.mengxiang@admerahealth.com</dev/null'.format(sample_sheet_name,flowcell_ID,sample_sheet_name,sample_sheet_name)
            subprocess.check_call(send_report,shell=True)
            print (send_report)


        '''
        #delete the duplicate Undetermined Fastq
        def find_files(pattern,path,run_id):
            result=[]
            for root,dirs,files in os.walk(path):
                if os.path.basename(root).startswith(run_id):
                    for name in files:
                        #print (name)
                        if fnmatch.fnmatch(name,pattern):
                            result.append(os.path.join(root,name))
            return result

        #sample_sheets=glob.glob('./SampleSheet*csv')
        sample_sheet_dict={}
        #m=re.match('SampleSheet-(.*)-(\d+)nt-(\d+)nt.csv',item)
        m=re.match('(SampleSheet-.*)-(\d+nt-\d+nt.csv)',sample_sheet_first)
        run_id=m.group(1)
        for sample_sheet in sample_sheets:
            sample_sheet_f=open(sample_sheet)
            for line in sample_sheet_f:
                if not line.startswith('Lane'):
                    continue
                if line.startswith('Lane'):
                    break
            for line in sample_sheet_f:
                items=line.split(',')
                lane_number=items[0]
                customer_id=items[1].split('-')[0]
                if lane_number not in sample_sheet_dict.keys():
                    sample_sheet_dict[lane_number]=[customer_id]
                elif customer_id not in sample_sheet_dict[lane_number] :
                    sample_sheet_dict[lane_number].append(customer_id)

        fastqs_need_to_delete=[]
        for key,value in sample_sheet_dict.items():
            if len(value) > 1:
                fastqs_need_to_delete.append('Undetermined_S0_L00%s_R1_001.fastq.gz' %key)
                fastqs_need_to_delete.append('Undetermined_S0_L00%s_R2_001.fastq.gz' %key)
                print (fastqs_need_to_delete)

        for need_delete_file in fastqs_need_to_delete:
            start_path="/mnt/novaoutput2/DeliverTeam_Output/6-Dmx-Fastq/"
            found_files=find_files(need_delete_file,start_path,run_id)

            for file_path in found_files:
                print(file_path)
                os.remove(file_path)
        '''



###check the run status, try to start demutiplex if run cycles finished.
def run_schedule():
    global all_run_finished_dmx

    run_schedule_info='Run_schedule.csv'
    with open(run_schedule_info,'r') as csvfile:
        reader=csv.reader(csvfile)
        run_schedule_list=list(reader)
    column_values=[row[2] for row in run_schedule_list]
    column_values+=[row[3] for row in run_schedule_list]
    ##check if all the runs finished demultiplex
    if '' not in column_values:
        all_run_finished_dmx=True
        #print (all_run_finished_dmx)
    for index,line in enumerate(run_schedule_list):
        if index==0:
            continue
        is_predem_done=line[2]
        is_dem_done=line[3]
        run_folder=line[1]
        if is_predem_done=='':
            print ('try to start to do predem method')
            print (run_folder)
            os.chdir(run_folder)
            nova_predem()

            os.chdir('/mnt/novaoutput2/si_command')
            #此处应该再修改成重新读取整个文件，根据index修改本行，其余保持不变
            with open(run_schedule_info,'r') as csvfile:
                reader_predem=csv.reader(csvfile)
                run_schedule_list_predem=list(reader_predem)

            if error_happened==True:
                run_schedule_list_predem[index][2]='error'
            elif error_happened==None:
                run_schedule_list_predem[index][2]=''
            elif error_happened==False:
                run_schedule_list_predem[index][2]='done'

            with open(run_schedule_info,'w',newline='') as csvfile:
                writer=csv.writer(csvfile)
                writer.writerows(run_schedule_list_predem)

        elif is_dem_done=='':
            print ('try to strart do dem method')
            #go to the run folder
            print (run_folder)
            os.chdir(run_folder)
            nova_dem()


            os.chdir('/mnt/novaoutput2/si_command')
            #此处应该再修改成重新读取整个文件，根据index修改本行，其余保持不变
            with open(run_schedule_info,'r') as csvfile:
                reader_dem=csv.reader(csvfile)
                run_schedule_list_dem=list(reader_dem)

            if error_happened==True:
                run_schedule_list_dem[index][3]='error'
            elif error_happened==None:
                run_schedule_list_dem[index][3]=''
            elif error_happened==False:
                run_schedule_list_dem[index][3]='done'


            with open(run_schedule_info,'w',newline='') as csvfile:
                writer=csv.writer(csvfile)
                writer.writerows(run_schedule_list_dem)

def check_missing_bcl(number_of_cycles):
    tree=ET.parse('RunInfo.xml')
    root=tree.getroot()
    '''
    for child in root:
        print(child.tag,child.attrib)
    total_cycle_number=0
    for reads in root.iter('Read'):
        reads_dict=reads.attrib
        NumCycles=int(reads_dict['NumCycles'])
        total_cycle_number+=NumCycles
    print (total_cycle_number)
    '''

    for flowcell in root.iter('FlowcellLayout'):
        lane_number_total_number=int(flowcell.attrib['LaneCount'])

    CopyComplete=True
    for lane_number in range(1,lane_number_total_number+1):
        bcl_folder='./Data/Intensities/BaseCalls/L00%d' %lane_number
        for cycle_number in range(1,number_of_cycles+1):
            bcl_folder='./Data/Intensities/BaseCalls/L00%d/C%d.1/' %(lane_number,cycle_number)
            bcl_file_1='L00%d_1.cbcl' %lane_number
            bcl_file_2='L00%d_2.cbcl' %lane_number

            if os.path.isfile(os.path.join(bcl_folder,bcl_file_1)) and os.path.isfile(os.path.join(bcl_folder,bcl_file_2)):
                pass
            else:
                CopyComplete=False
    return CopyComplete

            #elif not os.path.isfile(os.path.join(bcl_folder,bcl_file_1)):
            #    print ('lane%d,%s is missing' %(cycle_number,bcl_file_1))
            #elif not os.path.isfile(os.path.join(bcl_folder,bcl_file_2)):
            #    print ('lane%d,%s is missing' %(cycle_number,bcl_file_2))




# In[4]:


import xml.etree.ElementTree as ET
import os
import re
import subprocess
import glob
import csv
#import schedule
import time
import fnmatch


global all_run_finished_dmx
all_run_finished_dmx=False


while True:
    print ('job is onging')
    if all_run_finished_dmx==True:
        print ('All runs finished demtiplex')
        break
    run_schedule()
    time.sleep(300)
