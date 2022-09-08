rule proc_asassn:
     input:
        "src/data/ASASSN-21js/light_curve_f9818a9a-2dfc-4c33-ad97-826ae7b78a33.csv"
     output:
        "src/data/obs_ASASSN-21js_ASASSN.ecsv"
     conda:
        "environment.yml"
     script:
        "src/scripts/convert_asassn-21js.py"
