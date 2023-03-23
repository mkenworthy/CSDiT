rule proc_asassn21js:
     input:
        "src/data/ASASSN-21js/light_curve_f9818a9a-2dfc-4c33-ad97-826ae7b78a33.csv"
     output:
        "src/data/obs_ASASSN-21js_ASASSN.ecsv"
     conda:
        "environment.yml"
     script:
        "src/scripts/convert_asassn-21js.py"

rule proc_asassn21nn:
     input:
        "src/data/ASASSN-21nn/light_curve_49ee9f4d-675b-4cbd-8042-f13e024b41ba.csv"
     output:
        "src/data/obs_ASASSN-21nn_ASASSN.ecsv"
     conda:
        "environment.yml"
     script:
        "src/scripts/convert_asassn-21nn.py"

rule proc_asassn21sa:
     input:
        "src/data/ASASSN-21sa/light_curve_0d7aa71d-615a-42a3-836b-0cf23615ec7b.csv"
     output:
        "src/data/obs_ASASSN-21sa_ASASSN.ecsv"
     conda:
        "environment.yml"
     script:
        "src/scripts/convert_asassn-21sa.py"

rule proc_J1706:
     input:
        "src/data/J1706/light_curve_e7c6cf12-82c8-44b6-90ea-c9ee3c09e1ba.csv"
     output:
        "src/data/obs_J1706_ASASSN.ecsv"
     conda:
        "environment.yml"
     script:
        "src/scripts/convert_J1706.py"

rule proc_J1816:
     input:
        "src/data/J1816/light_curve_b1c2ca34-b951-4bf6-94aa-443a8dfe053c.csv"
     output:
        "src/data/obs_J1816_ASASSN.ecsv"
     conda:
        "environment.yml"
     script:
        "src/scripts/convert_J1816.py"

rule proc_asassn23ao:
     input:
        "src/data/ASASSN-23ao/light_curve_8607bbd1-ae8c-4b13-aa6d-18391e2a3d6b.csv"
     output:
        "src/data/obs_ASASSN-23ao_ASASSN.ecsv"
     conda:
        "environment.yml"
     script:
        "src/scripts/convert_asassn-23ao.py"