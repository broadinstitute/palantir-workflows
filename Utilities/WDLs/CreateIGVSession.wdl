version 1.0

workflow CreateIGVSession {
    input {
        Array[String]? bams
        Array[String]? vcfs
        Array[String]? interval_files

        String reference

        String output_name = "igv_session"
    }

    call MakeIGVXML {
        input:
            bams=bams,
            vcfs=vcfs,
            interval_files=interval_files,
            reference=reference,
            output_name=output_name
    }


    output {
        File igv_session = MakeIGVXML.igv_session
    }
}

task MakeIGVXML {
    input {
        Array[String]? bams
        Array[String]? vcfs
        Array[String]? interval_files

        String reference

        String output_name = "igv_session"
    }

    Int disk_size = 50
    Int cpu = 2
    Int memory = 8

    command <<<
        set -xueo pipefail

        python << CODE
        import xml.etree.ElementTree as ET

        # Parse WDL array into Python lists
        bams_list = ["~{default="" sep="\" , \"" bams}"]
        vcfs_list = ["~{default="" sep="\" , \"" vcfs}"]
        interval_files_list = ["~{default="" sep="\" , \"" interval_files}"]

        # Build XML file using nodes in ElementTree, with 'Session' at root
        session = ET.Element('Session')

        # Check reference sequence based on Terra defaults; else use provided reference
        if "~{reference}" == "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta" or "~{reference}" == "hg38":
            session.set('genome', 'hg38')
        elif "~{reference}" == "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta" or "~{reference}" == "hg19":
            session.set('genome', 'hg19')
        else:
            session.set('genome', "~{reference}")    # Requires IGV v2.1+

        # Create resources object to hold all files
        resources = ET.SubElement(session, 'Resources')

        # Add all resources from BAMs

        for bam in bams_list:
            bam_name = bam.split('/')[-1]

            bam_panel = ET.SubElement(session, 'Panel')
            bam_panel.set('name', f'{bam_name}-BAMPanel')

            res = ET.SubElement(resources, 'Resource')
            res.set('path', bam)
            res.set('type', 'bam')

            cov_track = ET.SubElement(bam_panel, 'Track')
            cov_track.set('id', f'{bam}_coverage')
            cov_track.set('name', f'{bam_name} Coverage')

            track = ET.SubElement(bam_panel, 'Track')
            track.set('id', bam)
            track.set('name', bam_name)

        # Add all resources from VCFs and VCF track
        vcf_panel = ET.SubElement(session, 'Panel')
        vcf_panel.set('name', 'VCFPanel')
        for vcf in vcfs_list:
            res = ET.SubElement(resources, 'Resource')
            res.set('path', vcf)
            res.set('type', 'vcf')

            vcf_name = vcf.split('/')[-1]
            track = ET.SubElement(vcf_panel, 'Track')
            track.set('id', vcf)
            track.set('name', vcf_name)

        # Create tracks for intervals and reference
        feature_panel = ET.SubElement(session, 'Panel')
        feature_panel.set('name', 'FeaturePanel')

        # Add each interval_file to resources and lower track
        for interval_file in interval_files_list:
            res = ET.SubElement(resources, 'Resource')
            res.set('path', interval_file)
            interval_type = interval_file.split('.')[-1]    # Check if .bed or .interval_list
            res.set('type', interval_type)

            interval_name = interval_file.split('/')[-1]
            track = ET.SubElement(feature_panel, 'Track')
            track.set('id', interval_file)
            track.set('name', interval_name)

        # Make reference sequence visible
        ref_seq = ET.SubElement(feature_panel, 'Track')
        ref_seq.set('id', 'Reference sequence')
        ref_seq.set('name', 'Reference sequence')

        # Collect session into ElementTree and write to XML file with header
        xml = ET.ElementTree(session)
        xml.write('~{output_name}.xml', encoding='UTF-8', xml_declaration=True)

        CODE

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory + "GB"
    }

    output {
        File igv_session = "~{output_name}.xml"
    }
}