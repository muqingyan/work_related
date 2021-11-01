version 1.0

workflow varScanTest{
    input{
        Array [File] Nbam_file
	    Array [File] Tbam_file
	    File reference
	    File reference_index
        File target_bed
    }

    scatter (bam in zip(Nbam_file,Tbam_file)){
        call samtoolsPileup as pileN{
            input:
                bam_file = bam.left,
                reference = reference,
                target_region = target_bed


        }
        call samtoolsPileup as pileT{
            input:
                bam_file = bam.right,
                reference = reference,
                target_region = target_bed


        }
        call varScanSomatic as varsomatic{
            input:
                Npileup = pileN.out_pileup,
                Tpileup = pileT.out_pileup
        }
        call varScanProcessSomatic as varprocess{
            input:
                snpRaw = varsomatic.out_snp,
                indelRaw = varsomatic.out_indel
        }

    }

    output{
        Array[File] snpVarscans = varsomatic.out_snp
        Array[File] indelVarscans = varsomatic.out_indel
        Array[File] varscanvcfs = varprocess.snv_vcf
#        Array[File] snp_Somatic = varprocess.out_snp_Somatic
#        Array[File] snp_Somatic_hc = varprocess.out_snp_Somatic_hc
#        Array[File] snp_Germline = varprocess.out_snp_Germline
#        Array[File] snp_Germline_hc = varprocess.out_snp_Germline_hc
#        Array[File] snp_LOH = varprocess.out_snp_LOH
#        Array[File] snp_LOH_hc = varprocess.out_snp_LOH_hc
#        Array[File] indel_Somatic = varprocess.out_indel_Somatic
#        Array[File] indel_Somatic_hc = varprocess.out_indel_Somatic_hc
#        Array[File] indel_Germline = varprocess.out_indel_Germline
#        Array[File] indel_Germline_hc = varprocess.out_indel_Germline_hc
#        Array[File] indel_LOH = varprocess.out_indel_LOH
#        Array[File] indel_LOH_hc = varprocess.out_indel_LOH_hc
    }

}

task samtoolsPileup{
    input{
        File bam_file
        File reference
        File target_region
        String sample_name = basename(bam_file,".bam")
    }
    command{
        echo ${sample_name}
        /usr/bin/samtools-1.9/samtools mpileup -q1 \
        -f ${reference} -l ${target_region} -o ${sample_name}.pileup ${bam_file}

    }
    runtime{
        docker:"gcr.io/natera-pfrm-4713/varscan:v.2.4.2"
        dx_instance_type: "mem1_ssd1_v2_x16"
    }
    output{
        File out_pileup = "${sample_name}.pileup"
    }

}

task varScanSomatic{
    input{
        File Npileup
        File Tpileup
        String sample_name = basename(Tpileup,".pileup")
    }
    command{
        java -jar /opt/varscan/VarScan.v2.4.2.jar somatic \
        ${Npileup} ${Tpileup} ${sample_name}.varscan2 --output-vcf 1 --somatic-p-value 0.25 --min-var-freq 0.01
    }
    runtime{
        docker:"gcr.io/natera-pfrm-4713/varscan:v.2.4.2"
        dx_instance_type: "mem1_ssd1_v2_x16"
    }
    output{
        File out_snp = "${sample_name}.varscan2.snp.vcf"
        File out_indel = "${sample_name}.varscan2.indel.vcf"
    }
}

task varScanProcessSomatic{
    input{
        File snpRaw
        File indelRaw
        String sample_name = sub(basename(snpRaw,".vcf"),".varscan2.snp","")
    }
    command<<<
        mkdir /varscan
        cp /home/dnanexus/inputs/*/~{sample_name}*vcf /varscan/
        echo "1st ls -- before varscanprocess"
        ls /varscan
        java -jar /opt/varscan/VarScan.v2.4.2.jar processSomatic /varscan/~{sample_name}.varscan2.snp.vcf --min-tumor-freq 0.01 --p-value 0.25
        java -jar /opt/varscan/VarScan.v2.4.2.jar processSomatic /varscan/~{sample_name}.varscan2.indel.vcf --min-tumor-freq 0.01 --p-value 0.25
        echo "2nd ls -- after process"
        ls /varscan
        egrep "PASS|#" /varscan/~{sample_name}.varscan2.snp.Somatic.vcf > /varscan/~{sample_name}.varscan2.snp.Somatic.pass.vcf
        awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' /varscan/~{sample_name}.varscan2.snp.Somatic.pass.vcf > /varscan/~{sample_name}.varscan2.snp.Somatic.pass.chr.vcf

        mv /varscan/~{sample_name}.varscan2.snp.Somatic.pass.chr.vcf /varscan/~{sample_name}-TUMOR-FXG.varscan2.somatic.snv.vcf
        echo "after rename"
        ls /varscan
#        tar -zvcf /home/dnanexus/outputs.tar.gz outputfiles
        cp /varscan/~{sample_name}-TUMOR-FXG.varscan2.somatic.snv.vcf /home/dnanexus/
        echo "ls /home/dnanexus"
        ls /home/dnanexus/

    >>>
    runtime{
        docker:"gcr.io/natera-pfrm-4713/varscan:v.2.4.2"
        dx_instance_type: "mem1_ssd1_v2_x16"
    }
    output{
        File snv_vcf = "${sample_name}-TUMOR-FXG.varscan2.somatic.snv.vcf"
#        File out_snp_Somatic = "${snpRaw}.Somatic"
#        File out_snp_Somatic_hc = "${snpRaw}.Somatic.hc"
#        File out_snp_Germline = "${snpRaw}.Germline"
#        File out_snp_Germline_hc = "${snpRaw}.Germline.hc"
#        File out_snp_LOH = "${snpRaw}.LOH"
#        File out_snp_LOH_hc = "${snpRaw}.LOH.hc"
#        File out_indel_Somatic = "${indelRaw}.Somatic"
#        File out_indel_Somatic_hc = "${indelRaw}.Somatic.hc"
#        File out_indel_Germline = "${indelRaw}.Germline"
#        File out_indel_Germline_hc = "${indelRaw}.Germline.hc"
#        File out_indel_LOH = "${indelRaw}.LOH"
#        File out_indel_LOH_hc = "${indelRaw}.LOH.hc"
    }

}