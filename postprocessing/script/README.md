Date: 2024/12/09

Created a script to read from a named directory, e.g. `field_new_200m_P0.008_RE40000_10_15_rand2_Htheta0.503`, and create a config file in the processed folder. Run it by calling

`python level1.py --path='field_new_200m_P0.008_RE40000_10_15_rand2_Htheta0.503' --label='Case1'`

Moving forward, a better way to streamline thing would be to create the config file with all parameters at the time of running the simulation. This eliminates the need to use long names for folders. 

Step1: running the simulation. We can all the executable with the parameters defined in the config file. I only know how to do this in python.

```
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Training ANNs. Read yaml from a directory and save models to the same directory.')
    parser.add_argument('--path', type=str, help='Directory of config file.')
    parser.add_argument('--rand', type=int, help='Random seeding (for data split).')
    args = parser.parse_args()
    
    with open(args.path + 'config.json', 'r') as f:
        config = json.load(f)
```

Step2: postprocessing (level1).