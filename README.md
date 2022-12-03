# Lego-Brick Approach to Coding for Network Communication
This repository provides MATLAB implementation of the Lego-brick design of codes for network communication. More details about the coding schemes are available here: [https://arxiv.org/abs/2211.07208](https://arxiv.org/abs/2211.07208)

## Overview of what is provided
Implementations of coding schemes are provided for the problems of:
- Lossy source coding
- Gelfand-Pinsker coding
- Marton coding for broadcast channels
- Coding for cloud radio access networks (downlink and uplink)

Polar codes with successive cancellation decoding are used as the constituent point-to-point channel codes. The rates of the polar codes are chosen "close" to their theoretical limits, as described in the paper. The key implementations of encoders/decoders that are provided in the src directory are re-used for multiple problems. Scripts that are specific to a coding problem are given in the corresponding directory. 

## Support
Any questions and/or suggestions from members of the information, coding and communication theory communities are warmly welcomed. Email: [nghaddar@ucsd.edu](mailto:nghaddar@ucsd.edu)

## BibTex Citation
If you use our work in a publication, we would appreciate using the following citations:

```
@misc{Ghaddar2022,
  url = {https://arxiv.org/abs/2211.07208},
  author = {Ghaddar, Nadim and Ganguly, Shouvik and Wang, Lele and Kim, Young-Han},  
  title = {A {Lego}-Brick Approach to Coding for Network Communication},  
  publisher = {arXiv},  
  year = {2022}
}

@misc{Ghaddar2022Github,
  author = {Ghaddar, Nadim and Ganguly, Shouvik and Wang, Lele and Kim, Young-Han},
  title = {{Lego}-Brick Approach to Coding for Network Communication},
  url = {https://github.com/nadimgh/lego-brick},
  version = {1.0.0},
  year = {2022}
}
```
