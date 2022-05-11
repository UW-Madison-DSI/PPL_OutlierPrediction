<!-- PROJECT LOGO -->

<br />
<p align="center">
  <div style="width:100%; text-align:center">
  <img src="./images/logo.svg" alt="Logo">
  </div>

  <h3 align="center">WasteWater  Outlier Detection Investigation Using PPL</h3>

  <p align="center">
    This project investigated the use of probabilistic programming techniques to detect outliers in Covid-19 wastewater data. 
  </p>
</p>

<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
  * [Built With](#built-with)
  * [Project Steps](#project-steps)
* [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
* [Usage](#usage)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)

<!-- ABOUT THE PROJECT -->
## About The Project

This project uses the Julia programming language and the Gen probabilistic modeling and inferencing framework to model the data and to identify outliers among the sample data points. 

<div class="figure" style="text-align:center">
    <img src="images/quadratic-spline.png" style="max-width:600px" />
    <div class="caption" style="font-style:italic">Sample Results</div>
</div>

### Project Steps
For ease of understanding, this project is broken up into a set of steps:

1. [Importing Data](./src/step01-importing-data/README.md)

In this step, we import our data from a file using Julia utilities such as dataframes. 

2. [Linear Model](./src/step02-linear-model/README.md)

Next, we use a linear probabilistic model to fit a line through the data points along with identifying outliers.

3. [Linear Spline](./src/step03-linear-spline/README.md)

Rather than trying to approximate our data with a single line, we break the data into groups and fit a set of linear segments through the sequence of groups.

4. [Linear Log Spline](./src/step04-linear-log-spline/README.md)

Since our data points encompass a large range of values and because epidemic trends frequently follow exponential trends, we transform the values to a log scale and fit linear segments to the data in the logarithmic space.

5. [Quadratic Spline](./src/step05-quadratic-spline/README.md)

In order to overcome the limitations of using a linear model to fit varying data, we use a quadratic model instead, fitting the data to a series of quadratic spline segments.

### Built With

* [Julia - The Julia Programming Language](https://julialang.org)
* [Gen - PPL Framework for Generative Modeling and Probabilistic Inference](https://www.gen.dev)
* [Jupyter - Notebooks for Interactive Computing](https://jupyter.org)

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Prerequisites

* Julia

To install Julia, download and run the appropriate installer:
[Julia Downloads](https://julialang.org/downloads/)

* Jupyter

To run Jupyter notebooks, follow the instructions provided:
[Jupyter Installation](https://jupyter.org/install)

### Installation

1. Clone the repo
```sh
git clone https://github.com/github_username/repo_name.git
```

<!-- USAGE EXAMPLES -->
## Usage

The project examples are provided in two different forms (1) as Julia script files that may be exectued  from the command line and (2) as Jupyter notebook files that may be viewed in a Jupyter editor.

1. Go to the desired project step
```sh
cd src/step01-importing-data
```


1. Run the Julia code from the command line OR open the notebook file in a Jupyter editor.
```sh
julia step01.jl
``` 
```sh
jupyter notebook step01.ipynb
```

<!-- CONTRIBUTING -->
## Contributing

We are not accepting contributions to this project at this time.

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.

<!-- CONTACT -->
## Contact

Steve Goldstein - (mailto:sgoldstein@wisc.edu) - email

Project Link: [https://github.com/github_username/repo_name](https://github.com/github_username/repo_name)

<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

* [Data Science Institute @ University of Wisconsin-Madison](http://datascience.wisc.edu)



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo.svg?style=flat-square
[contributors-url]: https://github.com/github_username/repo/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo.svg?style=flat-square
[forks-url]: https://github.com/github_username/repo/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo.svg?style=flat-square
[stars-url]: https://github.com/github_username/repo/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo.svg?style=flat-square
[issues-url]: https://github.com/github_username/repo/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo.svg?style=flat-square
[license-url]: https://github.com/github_username/repo/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=flat-square&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/github_username
[product-screenshot]: images/screenshot.png
