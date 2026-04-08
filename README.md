<!-- Title -->
<h1 align="left">
NeuralCrop.jl
</h1>

<!-- description -->
<p align="center">
  <strong> 🧑‍🌾 💧 ☀️ 🌾 A fast and flexible Julia framework for crop modelling, supporting fully process-based and hybrid simulations on CPUs and GPUs. </strong>
</p>

[![Build Status](https://github.com/yunan-l/NeuralCrop.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yunan-l/NeuralCrop.jl/actions/workflows/CI.yml?query=branch%3Amain)

NeuralCrop is a differentiable hybrid global gridded crop model (GGCM) that combines the strengths of the state-of-the-art GGCM [LPJmL](https://doi.org/10.5194/gmd-11-1343-2018) with machine learning approaches. By implementing components in a differentiable form to achieve seamless integration with machine learning methods, NeuralCrop enables end-to-end 'online training', with ML components optimized in tandem with the model dynamics.


## Citing

If you use NeuralCrop.jl in research or other activities, we would be grateful if you could mention NeuralCrop.jl and cite our paper:

> Lin, Yunan, et al. "NeuralCrop: Combining physics and machine learning for improved crop yield projections." arXiv preprint arXiv:2512.20177 (2025).

The bibtex entry for the paper is:

```bibtex
@article{lin2025neuralcrop,
  title={NeuralCrop: Combining physics and machine learning for improved crop yield projections},
  author={Lin, Yunan and Bathiany, Sebastian and Badri, Maha and Gelbrecht, Maximilian and Hess, Philipp and Groenke, Brian and Heinke, Jens and M{\"u}ller, Christoph and Boers, Niklas},
  journal={arXiv preprint arXiv:2512.20177},
  year={2025}
}
```