lefse_format_input.py lefse_input_anal_cat.txt lefse_input_anal_cat.in -c 2 -s -1 -u 1 -o 1000000
lefse_format_input.py lefse_input_anal_dog.txt lefse_input_anal_dog.in -c 2 -s -1 -u 1 -o 1000000
lefse_format_input.py lefse_input_throat_dog.txt lefse_input_throat_dog.in -c 2 -s -1 -u 1 -o 1000000
lefse_format_input.py lefse_input_throat_cat.txt lefse_input_throat_cat.in -c 2 -s -1 -u 1 -o 1000000

lefse_run.py lefse_input_anal_cat.in lefse_input_anal_cat.res
lefse_run.py lefse_input_anal_dog.in lefse_input_anal_dog.res
lefse_run.py lefse_input_throat_dog.in lefse_input_throat_dog.res
lefse_run.py lefse_input_throat_cat.in lefse_input_throat_cat.res

lefse_plot_res.py lefse_input_anal_cat.res lefse_input_anal_cat.png --dpi 300
lefse_plot_res.py lefse_input_anal_dog.res lefse_input_anal_dog.png --dpi 300
lefse_plot_res.py lefse_input_throat_dog.res lefse_input_throat_dog.png --dpi 300
lefse_plot_res.py lefse_input_throat_cat.res lefse_input_throat_cat.png --dpi 300
