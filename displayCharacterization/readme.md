The env.yml is what was used to create the environment in miniconda
The requirements.txt was then generated from that environment with `conda list --export > requirements.txt`

So, env.yml will be more flexible, and requirements.txt will be more precise (in terms of version numbers and dependencies).