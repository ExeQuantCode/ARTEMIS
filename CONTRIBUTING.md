# Contributing Guidelines

Thank you for considering contributing to the ARTEMIS project! We appreciate your time and effort. To ensure a smooth collaboration, please follow the guidelines provided below.

Please first discuss potential changes you wish to make to the project via issue (preferably), or email.

> And if you like the project, but just don't have time to contribute, that's fine. There are other easy ways to support the project and show your appreciation, which we would also be very happy about:
> - Star the project
> - Mention it on social media platforms
> - Refer this project in your project's readme
> - Mention the project at local meetups and tell your friends/colleagues

<!-- omit in toc -->

## Table of Contents
- [Code of Conduct](#code-of-conduct)
- [I Have a Question](#i-have-a-question)
- [I Want to Contribute](#i-want-to-contribute)
- [Reporting Bugs](#reporting-bugs)
- [Suggesting Enhancements](#suggesting-enhancements)
- [Code Style](#code-style)
- [Testing](#testing)
- [Documentation](#documentation)
- [Contact](#contact)
- [License](#license)


## Code of Conduct

This project and everyone participating in it is governed by the
[ARTEMIS Code of Conduct](CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code.
Please report unacceptable behavior to the [ARTEMIS developers](mailto:support@artemis-materials.co.uk?subject=ARTEMIS%20-%behaviour).

## I Have a Question

> If you want to ask a question, we assume that you have read the available [Documentation](README.md).

Before you ask a question, it is best to search for existing [Issues](https://github.com/ExeQuantCode/ARTEMIS/issues) that might help you.
In case you have found a suitable issue and still need clarification, you can write your question in this issue. It is also advisable to search the internet for answers first.

If you then still feel the need to ask a question and need clarification, we recommend the following:

- Open an [Issue](https://github.com/ExeQuantCode/ARTEMIS/issues/new).
- Provide as much context as you can about what you're running into.
- Provide project and platform versions (python, fortran, pip), depending on what seems relevant.

We will then take care of the issue as soon as possible.

## I Want To Contribute

> ### Legal Notice <!-- omit in toc -->
> When contributing to this project, you must agree that you have authored 100% of the content, that you have the necessary rights to the content and that the content you contribute may be provided under the project license.

### Reporting Bugs
If you encounter any issues or have suggestions for improvement, please open an [Issue](https://github.com/ExeQuantCode/ARTEMIS/issues/new) on the repository's issue tracker.

When reporting, please provide as much context as possible and describe the reproduction steps that someone else can follow to recreate the issue on their own.
This usually includes your code.
For good bug reports you should isolate the problem and create a reduced test case.



### Suggesting Enhancements

This section guides you through submitting an enhancement suggestion for ARTEMIS, **including completely new features and minor improvements to existing functionality**.
Following these guidelines will help maintainers and the community to understand your suggestion and find related suggestions.

<!-- omit in toc -->
#### Before Submitting an Enhancement

- Make sure that you are using the latest version.
- Read the compilable documentation carefully and find out if the functionality is already covered.
- Perform a [search](https://github.com/ExeQuantCode/ARTEMIS/issues) to see if the enhancement has already been suggested. If it has, add a comment to the existing issue instead of opening a new one.
- Find out whether your idea fits with the scope and aims of the project. It's up to you to make a strong case to convince the project's developers of the merits of this feature. Keep in mind that we want features that will be useful to the majority of our users and not just a small subset. If you're just targeting a minority of users, consider writing an add-on/plugin library.


### Contributing Code

This guide provides the recommended route to contributing to ARTEMIS:

1. Fork the repository.
2. Clone the forked repository to your local machine.
3. Create a new branch for your changes.
4. Make your changes and commit them.
5. Push the changes to your forked repository.
6. Open a pull request to the main repository.

When submitting your contributions, please ensure the following:
- Provide a clear and descriptive title for your pull request.
- Include a detailed description of the changes made.
- Reference any related issues or pull requests, if applicable.
- Write unit tests for your contributions
- Ensure all existing tests pass before submitting your changes.
- Update the documentation to reflect your changes, if necessary (i.e. through FORD style commenting).
- Provide examples and usage instructions, if applicable.

Follow the [Code Style](#code-style) when contributing code to this project to ensure compatibility and a uniform format to the project.


### Code Style
- Follow the code style and conventions set out in the [RAFFLE codebase](https://github.com/ExeQuantCode/RAFFLE). Moving forward, the [ARTEMIS codebase](https://github.com/ExeQuantCode/ARTEMIS) will be adopting this format and, as such, will transition all old commenting, documentation, and general code style to this form (this is likely to be a slow process). 
- Use meaningful variable and function names.
- Write clear and concise comments. For the Fortran library, use comments compatible with the [FORD Fortran Documenter](https://forddocs.readthedocs.io/en/stable/). The Fortran library does not yet support the FORD documenter, but the plan is to align it with the RAFFLE codebase and maintain the same level of documenter support. For the Python wrapper, use comments compatible with [pandoc](https://pandoc.org).



## Contact
If you have any questions or need further assistance, feel free to contact the [ARTEMIS developers](mailto:support@artemis-materials.co.uk?subject=ARTEMIS%20-%query).

## License
This project is licensed under the [GPL-3.0 License](LICENSE).

<!-- omit in toc -->
## Attribution
This guide is based on the **contributing-gen** and has been copied from the [graphstruc](https://github.com/nedtaylor/graphstruc) repository, with permission from the creator (Ned Taylor).
[Make your own](https://github.com/bttger/contributing-gen)!
