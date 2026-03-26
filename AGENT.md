# Agent Organization System

## Directory Structure

The project uses a three-tier organization system in the `agent/` directory:

- `agent/plans/` - High-level concepts and goals
- `agent/tasks/` - Actionable implementation items  
- `agent/docs/` - Current codebase reality documentation

## Naming Conventions

### Plans and Tasks
- Format: `NNN_descriptive_name_max_seven_words.md`
- Use auto-incrementing numbers (001, 002, etc.)
- Maximum 7 words in description
- Use underscores between words

### Docs
- Format: `subsystem_context.md`
- Named by subsystem and specific context they document
- Designed for easy reference by future agent sessions

## Content Structure

### Plans
- High-level concepts and goals
- Reference target docs for changes
- Never reference tasks (tasks reference plans)

### Tasks
- **Change-request**: Intent and high-level description
- **Plan reference**: Link to source plan (if applicable)
- **Steps**: Precise file and language structure changes needed
- **Doc references**: Target docs that will need updates

### Docs
- Document current codebase reality
- Reference source files only (never plans/tasks)
- Single source of truth about current state

## Workflow

1. **Planning** → Create plan in `agent/plans/`
2. **Task Generation** → Break plan into tasks in `agent/tasks/`
3. **Implementation** → Execute task steps
4. **Testing & Verification** → Confirm changes work
5. **Documentation Update** → Update relevant docs in `agent/docs/`

## Reference Flow

- Tasks MAY reference Plans (optional)
- Tasks/Plans MUST reference target Docs
- Docs MUST reference source files only
- No circular references allowed

This system ensures clear traceability from goals to implementation while maintaining docs as authoritative current state.
