# 2018-12-18 14:48:56

Had to create an instance at AWS Free Tier so students could play with bash. So,
here were the steps:

```bash
ssh -i ~/mykeypair.pem.txt ec2-user@ec2-3-16-151-211.us-east-2.compute.amazonaws.com
for i in {10..99}; do
    echo iamuser${i}:@@@@@:10${i}:513:Student Account ${i}:/home/iamuser${i}:/bin/bash >> batch_user_add.txt;
done
```

Then I had to set /etc/ssh/sshd_config to use Password authentication, and now
it works fine. I'll put the test files in the /tmp, and we can start playing
from there.



