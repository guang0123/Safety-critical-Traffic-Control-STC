# Safety-critical-Traffic-Control-STC

This provides the simulation code for Safety-critical Traffic Control proposed in
 [Safety-critical traffic control by connected automated vehicles](https://www.sciencedirect.com/science/article/abs/pii/S0968090X2300219X).

## Safe policy

We provides three commonly-used safe policy: constant time headway, time to collision, and stopping distance headway. Set `safepolicy.name` to choose different safe policies.

## Nominal controller

The designed STC synthesizes a safety-critical controller by minimizing deviation from a nominal controller. We use a feedback controller as one example nominal controller. Other pre-designed controllers can be used by re-writing the function `nominal_controller`.

## Citation

    @article{zhao2023safety,
      title={Safety-critical traffic control by connected automated vehicles},
      author={Zhao, Chenguang and Yu, Huan and Molnar, Tamas G},
      journal={Transportation research part C: emerging technologies},
      volume={154},
      pages={104230},
      year={2023},
      publisher={Elsevier}
    }
