# IMLCS 0.1.0-alpha

This is an experimental prototype for the research about
Incremental Multiple Longest Common Sub-Sequence ([IMLCS]). We include all
the propotypes we used for testing.

## Table of contents

- [Getting Started]
   - [Prerequisites]
   - [Installing]
   - [Running]
- [Contributing]
- [Versioning]
- [Authors]
- [License]
- [Acknowledgments]


## Getting Started

To get a copy of this software download or clone the GitHub repository.

Download:

```
wget https://github.com/LuisRusso-INESC-ID/IMLCS.git
```

Clone:

```
git clone git@github.com:LuisRusso-INESC-ID/IMLCS.git
```

### Prerequisites

This package was tested with [Debian], it will most likely work with other
Linux distributions and UNIX variants. Some version of the following
components must exist in the system.

* C compiler, [gcc] or [clang]
* [GNU Make]

### Installing

The makefile provides several compilation options, for debugging, profiling
etc.

To obtain a simple testing binary simply invoque the makefile in the `src`
directory.

```
cd src
make
```

If all went well this should produce a binary named `./project`. This is a
simple input shell. There are several altenative `main` functions to
execute unit tests.

The software in the sub-directories can be compiled in the same way. Be
careful that most of the files in there are actually links to files in the
`src` directory. In particular the `makefile` files. This makes it easy to
compile a consistent set of test binaries.

The algorithm in the `naive` directory is the classical dynamic programming
algorithm, similar classical [DP] algorithm but generalized for multiple
sequences. The algorithm in the `quickNaive` directory is an implementation
of the [QuickDP] algorithm, adapted for dynamic operations.

### Running

The input shell accepts lines that start by a letter, indicating the
respective command.

  - `X` terminates the program.

  - `K` followed by two numbers, defines the parameters of the input. The
    first number the number of strings and the second is the alphabet size.

  - `D` followed by a number, applies the `Pop()` operation to the string given
    by the number. The first string is number `0`.

  - `I` followed by a number and a letter, applies the `Append()` operation
    to the string given by the number.

The file `input` contains an example of a valid sequence of commands. This
is the example shown in the paper. The commands discussed in the paper are
the last two, before the `X` command.

You can test this input by executing:

```
./project < input
```

This outputs the number of `I` and `D` operations that were executed. The
binary has an internal `10` second timer. Defined by the constant
`TIME_LIMIT` in the makefile. The output is the number of operations
executed until the timer expires. In this case all the operations get
executed and the binary reaches the command `X`.

## Contributing

If you found this project useful please share it and the [IMLCS] article,
also you can create an [issue] with comments and suggestions, or email me to
[lmsrusso@gmail.com]

## Versioning

We use [SemVer] for versioning. For the versions available, see the [tags]
on this repository.

## Authors

* **Luís M. S. Russo** - *Initial work* - [LuisRusso]

See also the list of [contributors] who participated in this project.

## License

This project is licensed under the BSD 2-Clause "Simplified" License - see
the [LICENSE file] for details

## Acknowledgments

This prototype was develped for research that was supported by national
funds through Fundação para a Ciência e Tecnologia ([FCT]) through projects
[NGPHYLO] PTDC/CCI-BIO/29676/2017 and UID/CEC/50021/2019. Also Funded in part
by European Union’s Horizon 2020 research and innovation programme under
the Marie Sklodowska-Curie Actions grant agreement No 690941, [BIRDS].

* The participants of the [Lisbon String Masters 2018].
* Thanks to [PurpleBooth] for the [README-Template].
* The [grip] tool by [Joe Esposito] was very handy for producing this file.


[Getting Started]: #getting-started
[Prerequisites]: #prerequisites
[Installing]: #installing
[Running]: #running
[Contributing]: #contributing
[Versioning]: #versioning
[Authors]: #authors
[License]: #license
[Acknowledgments]: #acknowledgments

[IMLCS]: https://arxiv.org/abs/2005.02725
[Debian]: https://www.debian.org/
[gcc]: https://gcc.gnu.org/
[clang]: https://clang.llvm.org/
[GNU Make]: https://www.gnu.org/software/make/
[DP]: https://dl.acm.org/doi/10.1145/321796.321811
[QuickDP]: https://ieeexplore.ieee.org/document/5530316
[issue]: ../../issues
[lmsrusso@gmail.com]: mailto:lmsrusso@gmail.com
[SemVer]: http://semver.org/
[tags]: ../../tags
[LuisRusso]: https://github.com/LuisRusso
[contributors]: ../../contributors
[LICENSE file]: ./LICENSE
[FCT]: https://www.fct.pt/
[NGPHYLO]: https://thor.inesc-id.pt/ngphylo/
[BIRDS]: http://www.birdsproject.eu/
[README-Template]: https://gist.github.com/PurpleBooth/109311bb0361f32d87a2
[Lisbon String Masters 2018]: https://thor.inesc-id.pt/stringmasters/
[PurpleBooth]: https://gist.github.com/PurpleBooth
[grip]: https://github.com/joeyespo/grip
[Joe Esposito]: https://github.com/joeyespo
